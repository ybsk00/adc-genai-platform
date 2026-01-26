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
import os

from app.core.supabase import supabase

logger = logging.getLogger(__name__)

# ClinicalTrials.gov API v2 설정
API_BASE_URL = "https://clinicaltrials.gov/api/v2/studies"

# ADC 관련 검색어 (Split Mode - 안정성 및 진행률 표시 최적화)
ADC_SEARCH_TERMS = [
    'Antibody Drug Conjugate',
    '"Antibody-Drug Conjugate"',
    'ADC AND (Cancer OR Tumor OR Oncology OR Neoplasm)', # 노이즈 방지용 컨텍스트 추가
    'Immunoconjugate',
    'Trastuzumab Deruxtecan',
    'Enhertu',
    'DS-8201',
    'Sacituzumab Govitecan',
    'Trodelvy',
    'Brentuximab Vedotin',
    'Adcetris',
    'Ado-trastuzumab Emtansine',
    'Kadcyla',
    'T-DM1',
    'Polatuzumab Vedotin',
    'Polivy',
    'Loncastuximab Tesirine',
    'Zynlonta',
    'Tisotumab Vedotin',
    'Tivdak',
    'Mirvetuximab Soravtansine',
    'Elahere',
    'Datopotamab Deruxtecan',
    'Dato-DXd'
]

# 타겟 키워드 (drug name 추출용)
TARGET_KEYWORDS = ['adc', 'conjugate', 'antibody-drug', 'her2', 'trop2', 'bcma', 'dll3', 'nectin']


class BulkImporter:
    def __init__(self):
        self.total_imported = 0
        self.duplicates_skipped = 0
        self.errors = []

    def extract_study_info(self, study: dict) -> Dict[str, Any]:
        """임상시험 데이터에서 필요한 정보 추출 (ResultsSection 포함)"""
        protocol = study.get("protocolSection", {})
        results = study.get("resultsSection", {})
        
        id_module = protocol.get("identificationModule", {})
        status_module = protocol.get("statusModule", {})
        description_module = protocol.get("descriptionModule", {})
        arms_module = protocol.get("armsInterventionsModule", {})
        
        nct_id = id_module.get("nctId", "")
        title = id_module.get("officialTitle") or id_module.get("briefTitle", "No Title")
        
        # 약물 정보 추출
        interventions = arms_module.get("interventions", [])
        drug_names = [i.get("name", "") for i in interventions if i.get("type") in ["DRUG", "BIOLOGICAL"]]
        
        # --- 정량적 데이터 추출 (ResultsSection) ---
        outcome_data = self._parse_results_section(results)
        
        # 기본 정보 구성
        info = {
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
            "enrichment_source": "clinical_trials_api_v2",
            "confidence_score": 0.95 if outcome_data.get("has_results") else 0.5
        }
        
        # 정량 데이터 병합 (orr_pct, os_months, pfs_months, patient_count 등)
        info.update(outcome_data.get("metrics", {}))
        
        return info

    def _parse_results_section(self, results: dict) -> Dict[str, Any]:
        """ResultsSection에서 ORR, OS, PFS, Patient Count 등 추출"""
        metrics = {
            "orr_pct": None,
            "os_months": None,
            "pfs_months": None,
            "dor_months": None,
            "patient_count": None,
            "adverse_events_grade3_pct": None
        }
        
        if not results:
            return {"has_results": False, "metrics": metrics}
            
        # 1. 환자 수 (Participant Flow)
        participant_flow = results.get("participantFlowModule", {})
        groups = participant_flow.get("groups", [])
        if groups:
            try:
                # 'Total' 그룹이 명시적으로 있는지 확인 (중복 합산 방지)
                total_group = next((g for g in groups if "total" in g.get("title", "").lower()), None)
                if total_group:
                    metrics["patient_count"] = int(total_group.get("count", 0))
                else:
                    # Arms 합산
                    metrics["patient_count"] = sum(int(g.get("count", 0)) for g in groups if g.get("count"))
            except:
                pass

        # 2. 효능 지표 (Outcome Measures)
        outcome_measures = results.get("outcomeMeasuresModule", {}).get("outcomeMeasures", [])
        for measure in outcome_measures:
            title = measure.get("title", "").lower()
            unit = measure.get("unitOfMeasure", "").lower()
            classes = measure.get("classes", [])
            if not classes: continue
            
            # 첫 번째 클래스의 첫 번째 카테고리 데이터 사용 (단순화)
            categories = classes[0].get("categories", [])
            if not categories: continue
            
            measurements = categories[0].get("measurements", [])
            if not measurements: continue
            
            value_str = measurements[0].get("value")
            if not value_str: continue
            
            try:
                value = float(value_str)
                
                # 키워드 및 단위 검증 (AI 2차 검증 로직)
                if any(kw in title for kw in ["objective response rate", "orr", "overall response rate"]):
                    if self._verify_unit(unit, "percentage"):
                        metrics["orr_pct"] = value
                elif any(kw in title for kw in ["overall survival", "os"]) and "month" in title:
                    if self._verify_unit(unit, "months"):
                        metrics["os_months"] = value
                elif any(kw in title for kw in ["progression-free survival", "pfs"]) and "month" in title:
                    if self._verify_unit(unit, "months"):
                        metrics["pfs_months"] = value
                elif any(kw in title for kw in ["duration of response", "dor"]) and "month" in title:
                    if self._verify_unit(unit, "months"):
                        metrics["dor_months"] = value
            except ValueError:
                continue

        return {"has_results": True, "metrics": metrics}

    def _verify_unit(self, unit_str: str, expected_type: str) -> bool:
        """단위 검증: percentage vs months"""
        if not unit_str:
            return True # 단위가 없으면 일단 허용 (보수적 접근)
            
        unit_str = unit_str.lower()
        if expected_type == "percentage":
            # % 단위 확인 (percent, %, percentage)
            return any(u in unit_str for u in ["%", "percent"])
        elif expected_type == "months":
            # 시간 단위 확인 (month, week, year, day)
            return any(u in unit_str for u in ["month", "week", "year", "day"])
        return True

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
            return None 
        return None

    async def save_batch(self, batch: List[Dict[str, Any]], job_id: Optional[str] = None):
        """배치 단위로 DB에 저장 (Upsert 전략: Select -> Insert or Update)"""
        if not batch:
            return 0
        
        saved_count = 0
        for entry in batch:
            try:
                nct_id = entry.get("properties", {}).get("nct_id")
                if not nct_id:
                    continue
                
                # 1. 중복 체크 (nct_id 기준)
                existing = supabase.table("golden_set_library")\
                    .select("id")\
                    .eq("properties->>nct_id", nct_id)\
                    .execute()
                
                if existing.data:
                    # 2. 존재하면 업데이트 (Upsert 효과)
                    record_id = existing.data[0]['id']
                    # 업데이트할 필드만 선택 (기존 AI 분석 결과 등은 보존하고 싶을 수 있음)
                    # 여기서는 최신 데이터로 덮어쓰기 (properties 등)
                    update_data = {
                        "name": entry["name"],
                        "description": entry["description"],
                        "properties": entry["properties"],
                        "outcome_type": entry["outcome_type"],
                        # status나 ai_refined는 건드리지 않음 (이미 작업 중일 수 있으므로)
                    }
                    supabase.table("golden_set_library").update(update_data).eq("id", record_id).execute()
                    # logger.info(f"Updated existing record: {nct_id}")
                else:
                    # 3. 없으면 삽입
                    supabase.table("golden_set_library").insert(entry).execute()
                    saved_count += 1
                    self.total_imported += 1
                
            except Exception as e:
                logger.error(f"Save/Update error for {nct_id}: {e}")
                self.errors.append(str(e))
        
        return saved_count

    async def fetch_studies_generator(self, search_term: str, status_filter: List[str], page_size: int = 100, max_pages: int = 100, mode: str = "daily"):
        """API v2로 임상시험 데이터 조회 (제너레이터 방식 - 증분 처리 지원)"""
        next_page_token = None
        page = 0
        
        # 날짜 필터 설정 (Daily Sync용)
        last_update_date = None
        if mode == "daily":
            from datetime import timedelta
            yesterday = datetime.utcnow() - timedelta(days=1)
            last_update_date = yesterday.strftime("%Y-%m-%d")
            logger.info(f"📅 Daily Sync Mode: Fetching updates since {last_update_date}")
        
        # 프록시 설정
        proxy_url = None
        if os.getenv("PROXY_ENABLED", "").lower() == "true":
            host = os.getenv("PROXY_HOST")
            port = os.getenv("PROXY_PORT")
            user = os.getenv("PROXY_USERNAME")
            password = os.getenv("PROXY_PASSWORD")
            
            if host and port:
                if user and password:
                    proxy_url = f"http://{user}:{password}@{host}:{port}"
                else:
                    proxy_url = f"http://{host}:{port}"
                logger.info(f"🌐 Using proxy: {host}:{port}")

        async with aiohttp.ClientSession() as session:
            while page < max_pages:
                # API v2 syntax for date filtering
                query_term = search_term
                if last_update_date:
                    query_term = f"{search_term} AND AREA[LastUpdatePostDate]RANGE[{last_update_date},MAX]"
                
                params = {
                    "query.term": query_term,
                    "filter.overallStatus": ",".join(status_filter) if isinstance(status_filter, list) else status_filter,
                    "pageSize": page_size,
                    "format": "json"
                }
                
                if next_page_token:
                    params["pageToken"] = next_page_token
                
                # if last_update_date:
                #     params["filter.lastUpdatePostDate"] = last_update_date
                
                try:
                    request_kwargs = {
                        "params": params,
                        "timeout": aiohttp.ClientTimeout(total=120)
                    }
                    if proxy_url:
                        request_kwargs["proxy"] = proxy_url
                        
                    async with session.get(API_BASE_URL, **request_kwargs) as response:
                        if response.status != 200:
                            logger.error(f"API Error: HTTP {response.status}")
                            break
                        
                        data = await response.json()
                        studies = data.get("studies", [])
                        
                        if not studies:
                            break
                        
                        yield studies
                        
                        next_page_token = data.get("nextPageToken")
                        if not next_page_token:
                            break
                        
                        page += 1
                        await asyncio.sleep(0.5)  # Rate limiting
                        
                except Exception as e:
                    logger.error(f"Fetch error: {e}")
                    break

    async def run_import(self, job_id: Optional[str] = None, max_studies: int = 5000, mode: str = "daily"):
        """
        ClinicalTrials.gov API v2를 사용하여 ADC 임상시험 데이터 수집
        mode: 'daily' (기본값, 어제 이후 변경분) 또는 'full' (전체 덤프)
        """
        from app.api.scheduler import update_job_status, is_cancelled
        
        if job_id:
            await update_job_status(job_id, status="running")
        
        logger.info(f"🚀 [API v2 Importer] Starting ClinicalTrials.gov Split Search (Mode: {mode})...")
        
        try:
            # 상태 필터: 모든 상태 포함 (Broad Search)
            status_filter = [
                "COMPLETED", "TERMINATED", "WITHDRAWN", 
                "SUSPENDED", "RECRUITING", "ACTIVE_NOT_RECRUITING",
                "ENROLLING_BY_INVITATION", "NOT_YET_RECRUITING", "UNKNOWN"
            ]
            
            total_fetched = 0
            
            # Full 모드일 때는 페이지 제한을 넉넉하게
            # Daily 모드일 때는 적게 (어차피 날짜 필터로 걸러짐)
            max_pages_per_term = 100 if mode == "full" else 10
            
            for search_term in ADC_SEARCH_TERMS:
                if self.total_imported >= max_studies:
                    logger.info(f"✅ Reached max studies limit: {max_studies}")
                    break
                
                # 중단 요청 체크
                if job_id and await is_cancelled(job_id):
                    logger.info("Import cancelled by user")
                    await update_job_status(job_id, status="stopped")
                    return
                
                logger.info(f"🔍 Searching for: {search_term}")
                
                # 페이지 단위로 즉시 처리
                async for studies_page in self.fetch_studies_generator(
                    search_term=search_term,
                    status_filter=status_filter,
                    page_size=100,
                    max_pages=max_pages_per_term,
                    mode=mode
                ):
                    total_fetched += len(studies_page)
                    
                    batch = []
                    for study in studies_page:
                        if self.total_imported >= max_studies:
                            break
                        
                        entry = self.extract_study_info(study)
                        batch.append(entry)
                    
                    # 페이지 단위 즉시 저장
                    if batch:
                        saved = await self.save_batch(batch, job_id)
                        logger.info(f"💾 Page saved: {saved} new records (Term: {search_term[:20]}...)")
                        
                        # 진행률 실시간 업데이트
                        if job_id:
                            await update_job_status(
                                job_id, 
                                records_found=total_fetched,
                                records_drafted=self.total_imported
                            )
                    
                    if self.total_imported >= max_studies:
                        break
            
            logger.info(f"""
            ✅ [API v2 Import Complete]
            - Mode: {mode}
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
