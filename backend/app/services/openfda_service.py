import httpx
from typing import List, Dict, Any
from app.core.supabase import supabase
import logging

logger = logging.getLogger(__name__)

class OpenFDAService:
    BASE_URL = "https://api.fda.gov/drug/label.json"
    # AstraForge 2.0 Search Queries
    SEARCH_QUERIES = [
        'openfda.generic_name:"*vedotin"',
        'openfda.generic_name:"*deruxtecan"',
        'openfda.generic_name:"*govitecan"',
        'openfda.brand_name:"*ADC*"'
    ]

    async def fetch_approved_adcs(self, limit: int = 50) -> List[Dict[Any, Any]]:
        """승인된 ADC 라벨 정보 수집"""
        all_results = []
        async with httpx.AsyncClient(timeout=30.0) as client:
            for query in self.SEARCH_QUERIES:
                try:
                    params = {
                        "search": query,
                        "limit": limit
                    }
                    res = await client.get(self.BASE_URL, params=params)
                    if res.status_code == 200:
                        data = res.json()
                        results = data.get("results", [])
                        all_results.extend(results)
                        logger.info(f"FDA Query '{query}' found {len(results)} results.")
                except Exception as e:
                    logger.error(f"FDA Query Error '{query}': {e}")
        
        return all_results

    def extract_golden_info(self, label_data: Dict[Any, Any]) -> Dict[str, Any]:
        """FDA 라벨에서 golden_set_library 형식으로 정보 추출"""
        openfda = label_data.get("openfda", {})
        brand_name = openfda.get("brand_name", ["Unknown"])[0]
        generic_name = openfda.get("generic_name", ["Unknown"])[0]
        manufacturer = openfda.get("manufacturer_name", ["Unknown"])[0]
        
        # 독성 정보 (Boxed Warning) 추출
        boxed_warning = label_data.get("boxed_warning", ["No specific warning"])[0]
        
        return {
            "name": brand_name if brand_name != "Unknown" else generic_name,
            "category": "ADC",
            "description": f"FDA Approved ADC: {generic_name} by {manufacturer}.",
            "properties": {
                "generic_name": generic_name,
                "manufacturer": manufacturer,
                "boxed_warning": boxed_warning,
                "fda_label_id": label_data.get("id"),
                "approval_status": "Approved"
            },
            "outcome_type": "Success",
            "status": "approved",
            "enrichment_source": "openfda"
        }

    async def sync_to_db(self, job_id: Optional[str] = None):
        """FDA 데이터를 golden_set_library에 동기화"""
        from app.api.scheduler import sync_jobs
        
        labels = await self.fetch_approved_adcs()
        if job_id and job_id in sync_jobs:
            sync_jobs[job_id]["records_found"] = len(labels)
            
        drafted = 0
        for label in labels:
            try:
                # 중단 요청 체크
                if job_id and sync_jobs.get(job_id, {}).get("cancel_requested"):
                    if job_id in sync_jobs:
                        sync_jobs[job_id]["status"] = "stopped"
                    logger.info(f"Sync job {job_id} stopped by user.")
                    return

                golden_data = self.extract_golden_info(label)
                # 중복 체크 (generic_name 기준)
                existing = supabase.table("golden_set_library")\
                    .select("id")\
                    .eq("properties->>generic_name", golden_data["properties"]["generic_name"])\
                    .execute()
                
                if not existing.data:
                    supabase.table("golden_set_library").insert(golden_data).execute()
                    drafted += 1
                    if job_id and job_id in sync_jobs:
                        sync_jobs[job_id]["records_drafted"] = drafted
                    logger.info(f"✅ Synced Approved ADC: {golden_data['name']}")
            except Exception as e:
                logger.error(f"Sync Error for {label.get('id')}: {e}")
                if job_id and job_id in sync_jobs:
                    sync_jobs[job_id]["errors"].append(str(e))

openfda_service = OpenFDAService()
