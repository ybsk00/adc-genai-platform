import httpx
from datetime import datetime, timedelta
from typing import List, Dict, Any, Optional
from app.core.supabase import supabase
import logging
import asyncio

logger = logging.getLogger(__name__)

class OpenFDAService:
    BASE_URL = "https://api.fda.gov/drug/label.json"
    
    # ============================================================
    # EXPANDED SEARCH QUERIES - ADC 관련 모든 약물 포괄
    # ============================================================
    SEARCH_QUERIES = [
        # ADC 핵심 약물명 (Generic Name)
        'openfda.generic_name:vedotin',        # Brentuximab vedotin, Polatuzumab vedotin
        'openfda.generic_name:deruxtecan',     # Trastuzumab deruxtecan (Enhertu)
        'openfda.generic_name:govitecan',      # Sacituzumab govitecan (Trodelvy)
        'openfda.generic_name:emtansine',      # Trastuzumab emtansine (Kadcyla)
        'openfda.generic_name:ozogamicin',     # Gemtuzumab ozogamicin, Inotuzumab ozogamicin
        'openfda.generic_name:mafodotin',      # Belantamab mafodotin
        'openfda.generic_name:tesirine',       # Loncastuximab tesirine
        'openfda.generic_name:ravtansine',     # Mirvetuximab soravtansine
        'openfda.generic_name:tibsovo',
        
        # ADC 브랜드명 (Brand Name)
        'openfda.brand_name:ADCETRIS',         # Brentuximab vedotin
        'openfda.brand_name:KADCYLA',          # Trastuzumab emtansine
        'openfda.brand_name:ENHERTU',          # Trastuzumab deruxtecan
        'openfda.brand_name:TRODELVY',         # Sacituzumab govitecan
        'openfda.brand_name:PADCEV',           # Enfortumab vedotin
        'openfda.brand_name:BLENREP',          # Belantamab mafodotin
        'openfda.brand_name:ZYNLONTA',         # Loncastuximab tesirine
        'openfda.brand_name:ELAHERE',          # Mirvetuximab soravtansine
        'openfda.brand_name:MYLOTARG',         # Gemtuzumab ozogamicin
        'openfda.brand_name:BESPONSA',         # Inotuzumab ozogamicin
        'openfda.brand_name:POLIVY',           # Polatuzumab vedotin
        'openfda.brand_name:TIVDAK',           # Tisotumab vedotin
        
        # Antibody targets (항체 타겟)
        'openfda.generic_name:trastuzumab',    # HER2 targeting
        'openfda.generic_name:sacituzumab',    # Trop-2 targeting
        'openfda.generic_name:enfortumab',     # Nectin-4 targeting
        'openfda.generic_name:brentuximab',    # CD30 targeting
        'openfda.generic_name:polatuzumab',    # CD79b targeting
        'openfda.generic_name:tisotumab',      # Tissue Factor targeting
        'openfda.generic_name:gemtuzumab',     # CD33 targeting
        'openfda.generic_name:inotuzumab',     # CD22 targeting
        'openfda.generic_name:belantamab',     # BCMA targeting
        'openfda.generic_name:loncastuximab',  # CD19 targeting
        'openfda.generic_name:mirvetuximab',   # Folate receptor alpha
        
        # General conjugate terms
        'openfda.substance_name:"conjugate"',
        'openfda.pharm_class_epc:"antibody-drug conjugate"',
    ]

    async def fetch_all_approved_adcs(self, limit: int = 100, mode: str = "full") -> List[Dict[Any, Any]]:
        """
        OpenFDA에서 ADC 라벨 정보 수집 (Expanded)
        mode="full": 전체 데이터 수집 (Pagination)
        mode="daily": 최근 7일 업데이트된 데이터만 수집
        """
        all_results = []
        seen_ids = set()
        
        async with httpx.AsyncClient(timeout=60.0) as client:
            for idx, query in enumerate(self.SEARCH_QUERIES):
                skip = 0
                total_found = 0
                query_results = 0
                
                # Daily mode: Add date filter (7 days for better coverage)
                final_query = query
                if mode == "daily":
                    week_ago = (datetime.utcnow() - timedelta(days=7)).strftime("%Y%m%d")
                    final_query = f"{query}+AND+effective_time:[{week_ago}+TO+*]"

                # 1. Check Total Count first (Smart Check)
                try:
                    check_params = {"search": final_query, "limit": 1}
                    check_res = await client.get(self.BASE_URL, params=check_params)
                    if check_res.status_code == 200:
                        meta = check_res.json().get("meta", {})
                        total_found = meta.get("results", {}).get("total", 0)
                        logger.info(f"🔍 [{idx+1}/{len(self.SEARCH_QUERIES)}] Query '{query[:40]}...' found {total_found} records.")
                    elif check_res.status_code == 404:
                        # No results for this query
                        logger.debug(f"No results for query: {query[:40]}")
                        continue
                    else:
                        logger.warning(f"OpenFDA Check failed for '{query[:40]}': {check_res.status_code}")
                        continue
                except Exception as e:
                    logger.error(f"OpenFDA Check Error: {e}")
                    continue

                if total_found == 0:
                    continue

                # 2. Fetch Loop (All Pages)
                max_skip = 1000  # OpenFDA has a max skip limit
                while skip < min(total_found, max_skip):
                    try:
                        params = {
                            "search": final_query,
                            "limit": min(limit, 100),  # Max 100 per request
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
                                    query_results += 1
                            
                            if len(results) < limit:
                                break # End of results
                            
                            skip += limit
                            await asyncio.sleep(0.1) # Rate limiting
                        elif res.status_code == 404:
                            break
                        else:
                            logger.error(f"OpenFDA Fetch Error: {res.status_code}")
                            break
                    except Exception as e:
                        logger.error(f"OpenFDA Loop Exception: {e}")
                        break
                
                logger.info(f"   ↳ Added {query_results} unique records from this query.")
        
        logger.info(f"🎉 Total unique FDA labels fetched: {len(all_results)}")
        return all_results

    def extract_golden_info(self, label_data: Dict[Any, Any]) -> Dict[str, Any]:
        """FDA 라벨에서 golden_set_library 형식으로 정보 추출"""
        openfda = label_data.get("openfda", {})
        brand_name = openfda.get("brand_name", ["Unknown"])[0] if openfda.get("brand_name") else "Unknown"
        generic_name = openfda.get("generic_name", ["Unknown"])[0] if openfda.get("generic_name") else "Unknown"
        manufacturer = openfda.get("manufacturer_name", ["Unknown"])[0] if openfda.get("manufacturer_name") else "Unknown"
        
        # 독성 정보 (Boxed Warning) 추출
        boxed_warning = label_data.get("boxed_warning", ["No specific warning"])
        if isinstance(boxed_warning, list) and boxed_warning:
            boxed_warning = boxed_warning[0][:500]  # Truncate long warnings
        else:
            boxed_warning = "No specific warning"
        
        # Indication 추출
        indications = label_data.get("indications_and_usage", [""])
        indication_text = indications[0][:300] if isinstance(indications, list) and indications else ""
        
        return {
            "name": brand_name if brand_name != "Unknown" else generic_name,
            "category": "ADC",
            "description": f"FDA Approved: {generic_name} ({brand_name}) by {manufacturer}. {indication_text[:100]}",
            "properties": {
                "brand_name": brand_name,
                "generic_name": generic_name,
                "manufacturer": manufacturer,
                "boxed_warning": boxed_warning,
                "indication": indication_text,
                "fda_label_id": label_data.get("id"),
                "approval_status": "Approved",
                "last_updated": label_data.get("effective_time", datetime.utcnow().strftime("%Y%m%d"))
            },
            "outcome_type": "Success",
            "status": "draft",  # AI Refiner가 처리할 수 있도록 draft 상태
            "ai_refined": False,  # AI 정제 필요
            "enrichment_source": "open_fda_api"
        }

    async def sync_to_db(self, job_id: Optional[str] = None, mode: str = "full", limit: int = 100):
        """FDA 데이터를 golden_set_library에 동기화 (Smart Upsert)"""
        from app.api.scheduler import update_job_status, is_cancelled
        
        if job_id:
            await update_job_status(job_id, status="running")
            
        # Fetch data
        labels = await self.fetch_all_approved_adcs(limit=limit, mode=mode)
        
        if job_id:
            await update_job_status(job_id, records_found=len(labels))
        
        if not labels:
            logger.info("No FDA labels found to sync.")
            if job_id:
                await update_job_status(job_id, status="completed", message="No data found")
            return
            
        drafted = 0
        updated = 0
        errors = 0
        
        for idx, label in enumerate(labels):
            try:
                # Check cancellation
                if job_id and await is_cancelled(job_id):
                    await update_job_status(job_id, status="stopped")
                    logger.info(f"Sync job {job_id} stopped by user.")
                    return

                golden_data = self.extract_golden_info(label)
                drug_name = golden_data.get("name")
                
                if not drug_name or drug_name == "Unknown":
                    logger.warning(f"Skipping label with unknown name: {label.get('id')}")
                    continue
                
                # --- Smart Upsert Logic ---
                existing = supabase.table("golden_set_library")\
                    .select("id, properties")\
                    .eq("name", drug_name)\
                    .execute()
                
                if existing.data:
                    # Update existing record (Merge properties)
                    existing_record = existing.data[0]
                    existing_props = existing_record.get("properties", {}) or {}
                    new_props = golden_data.get("properties", {})
                    
                    # Merge: Keep existing data, add FDA data under 'fda_label' key
                    merged_props = existing_props.copy()
                    merged_props["fda_label"] = new_props
                    
                    supabase.table("golden_set_library")\
                        .update({
                            "properties": merged_props,
                            "updated_at": datetime.utcnow().isoformat(),
                        })\
                        .eq("id", existing_record["id"])\
                        .execute()
                    updated += 1
                    logger.debug(f"🔄 Updated: {drug_name}")
                else:
                    # Insert new record
                    try:
                        from app.services.chemical_resolver import chemical_resolver
                        smiles = chemical_resolver.fetch_verified_smiles(drug_name)
                        if smiles:
                            golden_data["smiles_code"] = smiles
                    except Exception as smiles_error:
                        logger.warning(f"SMILES lookup failed for {drug_name}: {smiles_error}")

                    supabase.table("golden_set_library").insert(golden_data).execute()
                    drafted += 1
                    logger.info(f"✅ Inserted: {drug_name}")

                # Update progress every 5 records
                if job_id and (drafted + updated) % 5 == 0:
                    await update_job_status(job_id, records_drafted=drafted + updated)
                    logger.info(f"📊 Progress: {drafted + updated}/{len(labels)} (Inserted: {drafted}, Updated: {updated})")

            except Exception as e:
                errors += 1
                logger.error(f"Sync Error for {label.get('id')}: {e}")

        logger.info(f"🎉 OpenFDA Sync Complete! Inserted: {drafted}, Updated: {updated}, Errors: {errors}")
        
        if job_id:
            await update_job_status(
                job_id, 
                status="completed", 
                records_drafted=drafted + updated,
                completed_at=datetime.utcnow().isoformat(), 
                message=f"Inserted: {drafted}, Updated: {updated}, Errors: {errors}"
            )

openfda_service = OpenFDAService()
