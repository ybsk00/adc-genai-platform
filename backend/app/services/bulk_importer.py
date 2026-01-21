"""
ClinicalTrials.gov API v2 Importer
덤프 파일 대신 공식 API v2를 사용하여 ADC 임상시험 데이터 수집
"""
import aiohttp
import asyncio
import logging
import json
from typing import Dict, Any, List, Optional
from datetime import datetime

from app.core.supabase import supabase

logger = logging.getLogger(__name__)

# ClinicalTrials.gov API v2 설정
API_BASE_URL = "https://clinicaltrials.gov/api/v2/studies"

# ADC 관련 검색어
ADC_SEARCH_TERMS = [
    "Antibody Drug Conjugate OR ADC OR Immunoconjugate", # 통합 쿼리 (Deep Scraping)
    "trastuzumab deruxtecan",
    "sacituzumab govitecan",
    "enfortumab vedotin",
    "HER2 conjugate",
    "TROP2 conjugate"
]

# 타겟 키워드 (drug name 추출용)
TARGET_KEYWORDS = ['adc', 'conjugate', 'antibody-drug', 'her2', 'trop2', 'bcma', 'dll3', 'nectin']


class BulkImporter:
    def __init__(self):
        self.total_imported = 0
        self.duplicates_skipped = 0
        self.errors = []

    def extract_study_info(self, study: dict) -> Dict[str, Any]:
        """임상시험 데이터에서 필요한 정보 추출"""
        protocol = study.get("protocolSection", {})
        id_module = protocol.get("identificationModule", {})
        status_module = protocol.get("statusModule", {})
        description_module = protocol.get("descriptionModule", {})
        arms_module = protocol.get("armsInterventionsModule", {})
        
        nct_id = id_module.get("nctId", "")
        title = id_module.get("officialTitle") or id_module.get("briefTitle", "No Title")
        
        # 약물 정보 추출
        interventions = arms_module.get("interventions", [])
        drug_names = [i.get("name", "") for i in interventions if i.get("type") == "DRUG"]
        
        return {
            "name": title[:200] if title else "Unknown",
            "category": "clinical_trial",
            "description": title,
            "properties": {
                "nct_id": nct_id,
                "phase": status_module.get("phases", []),
                "overall_status": status_module.get("overallStatus"),
                "why_stopped": status_module.get("whyStopped"),
                "brief_summary": description_module.get("briefSummary"),
                "drug_names": drug_names,
                "start_date": status_module.get("startDateStruct", {}).get("date"),
                "completion_date": status_module.get("completionDateStruct", {}).get("date"),
            },
            "status": "draft",
            "outcome_type": self._determine_outcome(status_module),
            "ai_refined": False,
            "enrichment_source": "clinical_trials_api_v2"
        }

    def _determine_outcome(self, status_module: dict) -> str:
        """상태를 기반으로 outcome_type 결정"""
        status = status_module.get("overallStatus", "").upper()
        why_stopped = status_module.get("whyStopped", "")
        
        if status == "COMPLETED":
            return "Success"
        elif status == "TERMINATED":
            if any(kw in why_stopped.lower() for kw in ["lack of efficacy", "futility", "not effective"]):
                return "Failure"
            return "Terminated"
        elif status in ["WITHDRAWN", "SUSPENDED"]:
            return "Terminated"
        elif status in ["RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION"]:
            # DB Constraint: outcome_type in ('Success', 'Failure', 'Terminated')
            # 'Ongoing' is not yet supported in DB constraint, so mapping to 'Unknown' (or null)
            # But 'Unknown' is not in the check list either?
            # Let's check the init.sql again: check (outcome_type in ('Success', 'Failure', 'Terminated'))
            # So 'Unknown' is ALSO not allowed if it's not null.
            # However, the column is nullable. So we should return None.
            return None 
        return None

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
                
                if existing.data:
                    self.duplicates_skipped += 1
                    continue
                
                # 새 레코드 저장
                supabase.table("golden_set_library").insert(entry).execute()
                saved_count += 1
                self.total_imported += 1
                
            except Exception as e:
                logger.error(f"Save error for {nct_id}: {e}")
                self.errors.append(str(e))
        
        return saved_count

    async def fetch_studies(self, search_term: str, status_filter: List[str], page_size: int = 100, max_pages: int = 100) -> List[dict]:
        """API v2로 임상시험 데이터 조회"""
        all_studies = []
        next_page_token = None
        page = 0
        
        async with aiohttp.ClientSession() as session:
            while page < max_pages:
                params = {
                    "query.term": search_term,
                    "filter.overallStatus": status_filter,
                    "pageSize": page_size,
                    "format": "json"
                }
                
                if next_page_token:
                    params["pageToken"] = next_page_token
                
                try:
                    async with session.get(API_BASE_URL, params=params, timeout=aiohttp.ClientTimeout(total=60)) as response:
                        if response.status != 200:
                            logger.error(f"API Error: HTTP {response.status}")
                            break
                        
                        data = await response.json()
                        studies = data.get("studies", [])
                        
                        if not studies:
                            break
                        
                        all_studies.extend(studies)
                        logger.info(f"📥 Fetched {len(studies)} studies (page {page + 1}, term: {search_term[:30]}...)")
                        
                        next_page_token = data.get("nextPageToken")
                        if not next_page_token:
                            break
                        
                        page += 1
                        await asyncio.sleep(0.5)  # Rate limiting
                        
                except Exception as e:
                    logger.error(f"Fetch error: {e}")
                    break
        
        return all_studies

    async def run_import(self, job_id: Optional[str] = None, max_studies: int = 500):
        """
        ClinicalTrials.gov API v2를 사용하여 ADC 임상시험 데이터 수집
        """
        from app.api.scheduler import update_job_status, is_cancelled
        
        if job_id:
            await update_job_status(job_id, status="running")
        
        logger.info("🚀 [API v2 Importer] Starting ClinicalTrials.gov data collection...")
        
        batch = []
        batch_size = 50
        
        try:
            # 상태 필터: 완료, 종료, 진행 중
            status_filters = [
                ["COMPLETED", "TERMINATED"],  # 결과가 있는 것
                ["ACTIVE_NOT_RECRUITING"],    # 진행 중이지만 데이터 있음
            ]
            
            total_fetched = 0
            
            for search_term in ADC_SEARCH_TERMS:
                if self.total_imported >= max_studies:
                    logger.info(f"✅ Reached max studies limit: {max_studies}")
                    break
                
                # 중단 요청 체크
                if job_id and await is_cancelled(job_id):
                    logger.info("Import cancelled by user")
                    await update_job_status(job_id, status="stopped")
                    return
                
                for status_filter in status_filters:
                    studies = await self.fetch_studies(
                        search_term=search_term,
                        status_filter=status_filter,
                        page_size=100,
                        max_pages=100 # Deep Scraping: 최대 100페이지(10,000건)까지 조회
                    )
                    
                    for study in studies:
                        if self.total_imported >= max_studies:
                            break
                        
                        entry = self.extract_study_info(study)
                        batch.append(entry)
                        total_fetched += 1
                        
                        # 배치 저장
                        if len(batch) >= batch_size:
                            saved = await self.save_batch(batch, job_id)
                            logger.info(f"💾 Batch saved: {saved} new records")
                            batch = []
                            
                            if job_id:
                                await update_job_status(
                                    job_id, 
                                    records_found=total_fetched,
                                    records_drafted=self.total_imported
                                )
            
            # 남은 배치 저장
            if batch:
                await self.save_batch(batch, job_id)
            
            logger.info(f"""
            ✅ [API v2 Import Complete]
            - Total Fetched: {total_fetched}
            - New Records: {self.total_imported}
            - Duplicates Skipped: {self.duplicates_skipped}
            - Errors: {len(self.errors)}
            """)
            
            if job_id:
                # 중복 스킵 정보를 errors 필드에 정보성으로 추가 (UI 표시용)
                completion_info = []
                if self.duplicates_skipped > 0:
                    completion_info.append(f"Info: {self.duplicates_skipped} records skipped (duplicate).")
                
                await update_job_status(
                    job_id,
                    status="completed",
                    records_found=total_fetched,
                    records_drafted=self.total_imported,
                    completed_at=datetime.utcnow().isoformat(),
                    errors=completion_info + self.errors # 기존 에러에 정보 추가
                )
                
        except Exception as e:
            error_msg = f"Import failed: {str(e)}"
            logger.error(error_msg)
            self.errors.append(error_msg)
            
            if job_id:
                await update_job_status(
                    job_id,
                    status="failed",
                    errors=self.errors
                )


bulk_importer = BulkImporter()
