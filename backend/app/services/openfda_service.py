import httpx
from datetime import datetime
from typing import List, Dict, Any, Optional
from app.core.supabase import supabase
import logging

logger = logging.getLogger(__name__)

class OpenFDAService:
    BASE_URL = "https://api.fda.gov/drug/label.json"
    # AstraForge 2.0 Search Queries (와일드카드 제거 - API 호환성)
    SEARCH_QUERIES = [
        'openfda.generic_name:vedotin',
        'openfda.generic_name:deruxtecan',
        'openfda.generic_name:govitecan',
        'openfda.generic_name:trastuzumab',
        'openfda.generic_name:sacituzumab',
        'openfda.generic_name:enfortumab'
    ]

    async def fetch_approved_adcs(self, limit: int = 50) -> List[Dict[Any, Any]]:
        """승인된 ADC 라벨 정보 수집"""
        all_results = []
        seen_ids = set()  # 중복 제거용
        
        async with httpx.AsyncClient(timeout=30.0) as client:
            for query in self.SEARCH_QUERIES:
                try:
                    params = {
                        "search": query,
                        "limit": limit
                    }
                    res = await client.get(self.BASE_URL, params=params)
                    logger.info(f"FDA Query '{query}' -> Status: {res.status_code}")
                    
                    if res.status_code == 200:
                        data = res.json()
                        results = data.get("results", [])
                        for r in results:
                            label_id = r.get("id")
                            if label_id and label_id not in seen_ids:
                                seen_ids.add(label_id)
                                all_results.append(r)
                        logger.info(f"FDA Query '{query}' found {len(results)} results (unique: {len(seen_ids)}).")
                    elif res.status_code == 404:
                        logger.info(f"FDA Query '{query}' returned 0 results.")
                    else:
                        logger.error(f"FDA Query '{query}' error: {res.status_code} - {res.text[:200]}")
                except Exception as e:
                    logger.error(f"FDA Query Exception '{query}': {e}")
        
        logger.info(f"Total unique FDA labels found: {len(all_results)}")
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
        """FDA 데이터를 golden_set_library에 동기화 (DB 기반)"""
        from app.api.scheduler import update_job_status, is_cancelled
        
        if job_id:
            await update_job_status(job_id, status="running")
            
        labels = await self.fetch_approved_adcs()
        if job_id:
            await update_job_status(job_id, records_found=len(labels))
            
        drafted = 0
        for label in labels:
            try:
                # 중단 요청 체크
                if job_id and await is_cancelled(job_id):
                    await update_job_status(job_id, status="stopped")
                    logger.info(f"Sync job {job_id} stopped by user.")
                    return

                golden_data = self.extract_golden_info(label)
                
                # SMILES 보강
                drug_name = golden_data.get("name")
                if drug_name:
                    from app.services.chemical_resolver import chemical_resolver
                    smiles = chemical_resolver.fetch_verified_smiles(drug_name)
                    if smiles:
                        golden_data["smiles_code"] = smiles

                # 중복 체크 (generic_name 기준)
                existing = supabase.table("golden_set_library")\
                    .select("id")\
                    .eq("properties->>generic_name", golden_data["properties"]["generic_name"])\
                    .execute()
                
                if not existing.data:
                    supabase.table("golden_set_library").insert(golden_data).execute()
                    drafted += 1
                    if job_id:
                        await update_job_status(job_id, records_drafted=drafted)
                    logger.info(f"✅ Synced Approved ADC: {golden_data['name']}")
            except Exception as e:
                logger.error(f"Sync Error for {label.get('id')}: {e}")
                if job_id:
                    # 에러 로그는 리스트로 관리 (실제 구현 시 DB 구조에 맞춰 조정)
                    pass

        if job_id:
            await update_job_status(job_id, status="completed", completed_at=datetime.utcnow().isoformat())

openfda_service = OpenFDAService()
