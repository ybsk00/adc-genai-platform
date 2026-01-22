import httpx
from datetime import datetime, timedelta
from typing import List, Dict, Any, Optional
from app.core.supabase import supabase
import logging
import asyncio

logger = logging.getLogger(__name__)

class OpenFDAService:
    BASE_URL = "https://api.fda.gov/drug/label.json"
    # AstraForge 2.0 Search Queries (Optimized with specific tags)
    # openfda.substance_name is more precise than general search
    SEARCH_QUERIES = [
        'openfda.substance_name:"antibody-drug conjugate"',
        'openfda.substance_name:"conjugate"',
        'openfda.generic_name:vedotin',
        'openfda.generic_name:deruxtecan',
        'openfda.generic_name:govitecan',
        'openfda.generic_name:trastuzumab',
        'openfda.generic_name:sacituzumab',
        'openfda.generic_name:enfortumab'
    ]

    async def fetch_all_approved_adcs(self, limit: int = 100, mode: str = "full") -> List[Dict[Any, Any]]:
        """
        OpenFDA에서 ADC 라벨 정보 수집
        mode="full": 전체 데이터 수집 (Pagination)
        mode="daily": 최근 업데이트된 데이터만 수집
        """
        all_results = []
        seen_ids = set()
        
        async with httpx.AsyncClient(timeout=60.0) as client:
            for query in self.SEARCH_QUERIES:
                skip = 0
                total_found = 0
                
                # Daily mode: Add date filter
                final_query = query
                if mode == "daily":
                    yesterday = (datetime.utcnow() - timedelta(days=1)).strftime("%Y-%m-%d")
                    final_query = f"{query} AND effective_time:[{yesterday} TO *]"

                # 1. Check Total Count first (Smart Check)
                try:
                    check_params = {"search": final_query, "limit": 1}
                    check_res = await client.get(self.BASE_URL, params=check_params)
                    if check_res.status_code == 200:
                        meta = check_res.json().get("meta", {})
                        total_found = meta.get("results", {}).get("total", 0)
                        logger.info(f"🔍 OpenFDA Query '{final_query}' found {total_found} total records.")
                        
                        # If total is huge, we might want to be careful, but for ADC it should be manageable.
                    else:
                        logger.warning(f"OpenFDA Check failed for '{final_query}': {check_res.status_code}")
                        continue
                except Exception as e:
                    logger.error(f"OpenFDA Check Error: {e}")
                    continue

                if total_found == 0:
                    continue

                # 2. Fetch Loop
                while True:
                    try:
                        params = {
                            "search": final_query,
                            "limit": limit,
                            "skip": skip
                        }
                        res = await client.get(self.BASE_URL, params=params)
                        
                        if res.status_code == 200:
                            data = res.json()
                            results = data.get("results", [])
                            if not results:
                                break
                                
                            for r in results:
                                label_id = r.get("id")
                                if label_id and label_id not in seen_ids:
                                    seen_ids.add(label_id)
                                    all_results.append(r)
                            
                            logger.info(f"Fetched {len(results)} items (Skip: {skip}) for '{final_query}'")
                            
                            if len(results) < limit:
                                break # End of results
                            
                            skip += limit
                            await asyncio.sleep(0.2) # Rate limiting
                        else:
                            logger.error(f"OpenFDA Fetch Error: {res.status_code} - {res.text[:200]}")
                            break
                    except Exception as e:
                        logger.error(f"OpenFDA Loop Exception: {e}")
                        break
        
        logger.info(f"Total unique FDA labels fetched: {len(all_results)}")
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
                "approval_status": "Approved",
                "last_updated": label_data.get("effective_time", datetime.utcnow().strftime("%Y%m%d"))
            },
            "outcome_type": "Success",
            "status": "approved",
            "enrichment_source": "open_fda_api"
        }

    async def sync_to_db(self, job_id: Optional[str] = None, mode: str = "daily", limit: int = 100):
        """FDA 데이터를 golden_set_library에 동기화 (Smart Upsert)"""
        from app.api.scheduler import update_job_status, is_cancelled
        
        if job_id:
            await update_job_status(job_id, status="running")
            
        # Fetch data
        labels = await self.fetch_all_approved_adcs(limit=limit, mode=mode)
        
        if job_id:
            await update_job_status(job_id, records_found=len(labels))
            
        drafted = 0
        updated = 0
        
        for label in labels:
            try:
                # Check cancellation
                if job_id and await is_cancelled(job_id):
                    await update_job_status(job_id, status="stopped")
                    logger.info(f"Sync job {job_id} stopped by user.")
                    return

                golden_data = self.extract_golden_info(label)
                drug_name = golden_data.get("name")
                
                # --- Smart Upsert Logic ---
                # 1. Try to find existing record by Name
                existing = supabase.table("golden_set_library")\
                    .select("*")\
                    .eq("name", drug_name)\
                    .execute()
                
                if existing.data:
                    # Update existing record (Merge properties)
                    existing_record = existing.data[0]
                    existing_props = existing_record.get("properties", {})
                    new_props = golden_data.get("properties", {})
                    
                    # Deep merge properties (Simple version: Update existing keys with new ones, keep others)
                    # For a true deep merge, we'd need a recursive function, but for now flat merge of top keys is usually enough
                    # or we specifically merge 'fda_data' into a separate key if we wanted isolation.
                    # User requested: "Clinical info and FDA info merged with different keys"
                    # Let's put FDA data into a 'fda_label' key inside properties to avoid overwriting clinical data
                    
                    merged_props = existing_props.copy()
                    merged_props["fda_label"] = new_props # Store FDA specific data in a sub-object
                    # Also update top-level fields if they are missing in existing
                    if not existing_props.get("manufacturer"):
                        merged_props["manufacturer"] = new_props.get("manufacturer")
                    
                    supabase.table("golden_set_library")\
                        .update({
                            "properties": merged_props,
                            "updated_at": datetime.utcnow().isoformat(),
                            # Don't overwrite enrichment_source if it's already set to something else (e.g. clinical_trials)
                            # But maybe we want to append? For now, leave it.
                        })\
                        .eq("id", existing_record["id"])\
                        .execute()
                    updated += 1
                    logger.info(f"🔄 Updated existing drug: {drug_name}")
                else:
                    # Insert new record
                    # SMILES enrichment for new records
                    if drug_name:
                        from app.services.chemical_resolver import chemical_resolver
                        smiles = chemical_resolver.fetch_verified_smiles(drug_name)
                        if smiles:
                            golden_data["smiles_code"] = smiles

                    supabase.table("golden_set_library").insert(golden_data).execute()
                    drafted += 1
                    logger.info(f"✅ Inserted new FDA drug: {drug_name}")

                if job_id and (drafted + updated) % 10 == 0:
                    await update_job_status(job_id, records_drafted=drafted) # We can track updated too if we add a field

            except Exception as e:
                logger.error(f"Sync Error for {label.get('id')}: {e}")
                # Continue to next item

        if job_id:
            await update_job_status(job_id, status="completed", completed_at=datetime.utcnow().isoformat(), message=f"Inserted: {drafted}, Updated: {updated}")

openfda_service = OpenFDAService()
