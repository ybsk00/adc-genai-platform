"""
ClinicalTrials.gov Bulk Importer Service
API 차단 우회를 위한 JSON 덤프 다운로드 및 일괄 적재
"""
import asyncio
import aiohttp
import zipfile
import io
import json
import logging
from typing import List, Dict, Any, Optional
from datetime import datetime

from app.core.supabase import supabase

logger = logging.getLogger(__name__)

# ClinicalTrials.gov 전체 데이터 덤프 URL
DUMP_URL = "https://clinicaltrials.gov/AllPublicJSON.zip"

# ADC 관련 키워드 필터
TARGET_KEYWORDS = [
    "antibody drug conjugate", "antibody-drug conjugate", "adc",
    "her2", "trop2", "egfr", "cd19", "cd22", "cd33", "cd30",
    "nectin-4", "bcma", "folate receptor",
    "vedotin", "deruxtecan", "govitecan", "emtansine", "ozogamicin",
    "mafodotin", "trastuzumab", "sacituzumab", "enfortumab",
    "polatuzumab", "brentuximab", "inotuzumab", "gemtuzumab"
]

class BulkImporter:
    def __init__(self):
        self.total_processed = 0
        self.total_imported = 0
        self.errors = []

    def is_adc_related(self, study_data: dict) -> bool:
        """ADC 관련 임상시험인지 필터링"""
        try:
            # 전체 JSON을 문자열로 변환하여 키워드 검색
            full_text = json.dumps(study_data).lower()
            return any(keyword in full_text for keyword in TARGET_KEYWORDS)
        except Exception:
            return False

    def extract_study_info(self, study_data: dict) -> Dict[str, Any]:
        """임상시험 데이터에서 필요한 정보 추출"""
        protocol = study_data.get("protocolSection", {})
        id_module = protocol.get("identificationModule", {})
        status_module = protocol.get("statusModule", {})
        description_module = protocol.get("descriptionModule", {})
        
        nct_id = id_module.get("nctId", "")
        title = id_module.get("officialTitle") or id_module.get("briefTitle", "No Title")
        
        return {
            "name": title[:200] if title else "Unknown",
            "category": "clinical_trial",
            "description": title,
            "properties": {
                "nct_id": nct_id,
                "phase": status_module.get("phase"),
                "overall_status": status_module.get("overallStatus"),
                "why_stopped": status_module.get("whyStopped"),
                "brief_summary": description_module.get("briefSummary"),
                "raw_data": study_data  # 전체 원본 데이터 저장 (AI Refiner용)
            },
            "status": "draft",
            "outcome_type": "Unknown",  # AI Refiner가 나중에 채움
            "ai_refined": False,  # 미정제 상태
            "enrichment_source": "clinical_trials_bulk"
        }

    async def save_batch(self, batch: List[Dict[str, Any]], job_id: Optional[str] = None):
        """배치 단위로 DB에 저장 (upsert)"""
        if not batch:
            return 0
        
        saved_count = 0
        for entry in batch:
            try:
                nct_id = entry.get("properties", {}).get("nct_id")
                if not nct_id:
                    continue
                
                # 중복 체크 (nct_id 기준)
                existing = supabase.table("golden_set_library")\
                    .select("id")\
                    .eq("properties->>nct_id", nct_id)\
                    .execute()
                
                if not existing.data:
                    supabase.table("golden_set_library").insert(entry).execute()
                    saved_count += 1
            except Exception as e:
                self.errors.append(f"Save error for {entry.get('properties', {}).get('nct_id')}: {str(e)[:100]}")
        
        return saved_count

    async def run_import(self, job_id: Optional[str] = None, max_studies: int = 5000):
        """
        Bulk Import 실행
        - JSON 덤프 다운로드 (스트리밍)
        - ADC 관련 데이터 필터링
        - 일괄 DB 적재
        """
        from app.api.scheduler import update_job_status, is_cancelled
        
        if job_id:
            await update_job_status(job_id, status="running")
        
        logger.info("🚀 [Bulk Importer] Starting ClinicalTrials.gov dump download...")
        
        batch = []
        batch_size = 100
        
        try:
            # ClinicalTrials.gov 403 방지를 위한 완전한 브라우저 헤더
            headers = {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
                'Accept': 'application/zip,application/octet-stream,*/*',
                'Accept-Language': 'en-US,en;q=0.9',
                'Accept-Encoding': 'gzip, deflate, br',
                'Connection': 'keep-alive',
                'Referer': 'https://clinicaltrials.gov/',
            }
            
            async with aiohttp.ClientSession(headers=headers) as session:
                logger.info(f"📥 Downloading from {DUMP_URL}...")
                
                async with session.get(DUMP_URL, timeout=aiohttp.ClientTimeout(total=3600)) as response:
                    if response.status != 200:
                        error_msg = f"Download failed: HTTP {response.status}"
                        logger.error(error_msg)
                        if job_id:
                            await update_job_status(job_id, status="failed", errors=[error_msg])
                        return
                    
                    # 스트리밍 다운로드 (메모리 효율적: 2GB → ~10MB)
                    import tempfile
                    import os
                    temp_path = None
                    
                    try:
                        logger.info("📦 Streaming ZIP file to temp... (this may take several minutes)")
                        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.zip')
                        temp_path = temp_file.name
                        downloaded_size = 0
                        
                        async for chunk in response.content.iter_chunked(1024 * 1024):  # 1MB chunks
                            temp_file.write(chunk)
                            downloaded_size += len(chunk)
                            if downloaded_size % (100 * 1024 * 1024) == 0:  # Log every 100MB
                                logger.info(f"📥 Downloaded {downloaded_size // (1024*1024)} MB...")
                        
                        temp_file.close()
                        logger.info(f"✅ Download complete: {downloaded_size // (1024*1024)} MB")
                        
                        if job_id:
                            await update_job_status(job_id, records_found=0)
                        
                        logger.info("📂 Extracting and parsing JSON files...")
                        
                        with zipfile.ZipFile(temp_path) as z:
                            json_files = [f for f in z.namelist() if f.endswith('.json')]
                            total_files = len(json_files)
                            logger.info(f"Found {total_files} JSON files in archive")
                        
                            for idx, filename in enumerate(json_files):
                                # 중단 요청 체크
                                if job_id and await is_cancelled(job_id):
                                    logger.info("Import cancelled by user")
                                    await update_job_status(job_id, status="stopped")
                                    return
                                
                                # 최대 수집 개수 체크
                                if self.total_imported >= max_studies:
                                    logger.info(f"Reached max studies limit: {max_studies}")
                                    break
                                
                                try:
                                    with z.open(filename) as f:
                                        study_data = json.load(f)
                                        self.total_processed += 1
                                        
                                        # ADC 관련 데이터만 필터링
                                        if self.is_adc_related(study_data):
                                            entry = self.extract_study_info(study_data)
                                            batch.append(entry)
                                            
                                            # 배치가 차면 저장
                                            if len(batch) >= batch_size:
                                                saved = await self.save_batch(batch, job_id)
                                                self.total_imported += saved
                                                batch = []
                                                
                                                if job_id:
                                                    await update_job_status(
                                                        job_id, 
                                                        records_found=self.total_processed,
                                                        records_drafted=self.total_imported
                                                    )
                                                
                                                logger.info(f"Progress: {self.total_processed}/{total_files} files, {self.total_imported} ADC studies imported")
                                
                                except json.JSONDecodeError:
                                    continue
                                except Exception as e:
                                    self.errors.append(str(e)[:100])
                        
                            # 남은 배치 저장
                            if batch:
                                saved = await self.save_batch(batch, job_id)
                                self.total_imported += saved
                    
                    finally:
                        # 임시 파일 반드시 삭제 (디스크 공간 확보)
                        if temp_path and os.path.exists(temp_path):
                            os.unlink(temp_path)
                            logger.info(f"🗑️ Temp file deleted: {temp_path}")
            
            # 완료
            logger.info(f"🎉 Import Complete! Total: {self.total_imported} ADC studies from {self.total_processed} files")
            
            if job_id:
                await update_job_status(
                    job_id,
                    status="completed",
                    records_found=self.total_processed,
                    records_drafted=self.total_imported,
                    completed_at=datetime.utcnow().isoformat(),
                    errors=self.errors[:20]
                )
        
        except Exception as e:
            logger.error(f"Bulk Import Error: {e}")
            if job_id:
                await update_job_status(job_id, status="failed", errors=[str(e)])

# 싱글톤 인스턴스
bulk_importer = BulkImporter()
