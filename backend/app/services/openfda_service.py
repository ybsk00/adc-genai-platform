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
    
    # ============================================================
    # BROADENED SEARCH QUERIES - Full-text 검색 (와일드카드 제거)
    # ============================================================
    BROADENED_QUERIES = [
        # 라벨 전문(Full Text) 검색 - ADC 관련 키워드
        'description:"antibody-drug conjugate"',
        'description:"antibody drug conjugate"',
        'description:"immunoconjugate"',
        'indications_and_usage:"conjugate"',
        'mechanism_of_action:"linker"',
        'mechanism_of_action:"payload"',
        'mechanism_of_action:"cytotoxic"',
        
        # 약리학적 분류(Pharmacologic Class) 검색
        'openfda.pharm_class_moa:"monoclonal antibody"',
        'openfda.pharm_class_cs:"antineoplastic"',
        
        # 암 면역치료제 관련 키워드
        'indications_and_usage:"cancer"',
        'indications_and_usage:"lymphoma"',
        'indications_and_usage:"leukemia"',
        'indications_and_usage:"myeloma"',
        'indications_and_usage:"carcinoma"',
    ]
    
    def _get_clinical_trial_drug_names(self) -> List[str]:
        """ClinicalTrials에서 이미 수집된 약물명 리스트 추출 (Targeted Search용)"""
        try:
            result = supabase.table("golden_set_library")\
                .select("name")\
                .eq("enrichment_source", "clinical_trials_api_v2")\
                .execute()
            
            if result.data:
                names = [r["name"] for r in result.data if r.get("name")]
                logger.info(f"📋 Found {len(names)} drug names from ClinicalTrials for targeted search")
                return names
            return []
        except Exception as e:
            logger.error(f"Failed to fetch ClinicalTrials drug names: {e}")
            return []
    
    def _generate_date_ranges(self, start_year: int = 2000) -> List[str]:
        """연도별 날짜 범위 생성 (Skip 1,000 한계 우회용)"""
        current_year = datetime.utcnow().year
        ranges = []
        for year in range(start_year, current_year + 1):
            start = f"{year}0101"
            end = f"{year}1231"
            ranges.append(f"effective_time:[{start}+TO+{end}]")
        return ranges
    
    # ============================================================
    # Rate Limiting & Retry Configuration
    # ============================================================
    REQUEST_DELAY = 1.0  # 요청 간 최소 지연 시간 (초)
    MAX_RETRIES = 3      # 최대 재시도 횟수
    INITIAL_BACKOFF = 2  # 첫 재시도 대기 시간 (초)
    
    async def _fetch_with_retry(self, client: httpx.AsyncClient, url: str, params: dict, max_retries: int = 3) -> Optional[httpx.Response]:
        """지수 백오프(Exponential Backoff)를 사용한 재시도 로직"""
        for attempt in range(max_retries):
            try:
                res = await client.get(url, params=params)
                
                # 성공 또는 404 (데이터 없음)
                if res.status_code in [200, 404]:
                    return res
                
                # 502, 503, 429 등 재시도 가능한 에러
                if res.status_code in [429, 500, 502, 503, 504]:
                    wait_time = self.INITIAL_BACKOFF * (2 ** attempt)  # 2, 4, 8초
                    logger.warning(f"⚠️ OpenFDA {res.status_code} error. Retry {attempt + 1}/{max_retries} in {wait_time}s...")
                    await asyncio.sleep(wait_time)
                    continue
                
                # 다른 에러는 즉시 반환
                return res
                
            except (httpx.TimeoutException, httpx.ConnectError) as e:
                wait_time = self.INITIAL_BACKOFF * (2 ** attempt)
                logger.warning(f"⚠️ Connection error: {e}. Retry {attempt + 1}/{max_retries} in {wait_time}s...")
                await asyncio.sleep(wait_time)
                continue
            except Exception as e:
                logger.error(f"❌ Unexpected error in fetch: {e}")
                return None
        
        logger.error(f"❌ Max retries ({max_retries}) exceeded.")
        return None

    async def fetch_all_approved_adcs(self, limit: int = 100, mode: str = "full", use_broad_search: bool = True, use_targeted_search: bool = True) -> List[Dict[Any, Any]]:
        """
        OpenFDA에서 ADC 라벨 정보 수집 (Rate Limited + Retry)
        mode="full": 전체 데이터 수집
        mode="daily": 최근 7일 업데이트된 데이터만 수집
        use_broad_search: True면 BROADENED_QUERIES도 함께 사용
        use_targeted_search: True면 ClinicalTrials 약물명으로 정밀 검색 추가
        """
        all_results = []
        seen_ids = set()
        
        # 1. 기본 쿼리 목록 합성
        queries_to_run = self.SEARCH_QUERIES.copy()
        if use_broad_search:
            queries_to_run.extend(self.BROADENED_QUERIES)
        
        # 2. Targeted Search: ClinicalTrials 약물명 기반 쿼리 추가
        if use_targeted_search:
            drug_names = self._get_clinical_trial_drug_names()
            for name in drug_names:
                # 약물명을 brand_name과 generic_name으로 모두 검색
                safe_name = name.replace('"', '\\"')  # 따옴표 이스케이프
                queries_to_run.append(f'openfda.brand_name:"{safe_name}"')
                queries_to_run.append(f'openfda.generic_name:"{safe_name}"')
        
        total_queries = len(queries_to_run)
        logger.info(f"🚀 Starting OpenFDA fetch with {total_queries} queries (mode={mode}, broad={use_broad_search}, targeted={use_targeted_search})")
        logger.info(f"⏱️ Rate limiting: {self.REQUEST_DELAY}s delay, {self.MAX_RETRIES} retries with exponential backoff")
        
        async with httpx.AsyncClient(timeout=60.0) as client:
            for idx, query in enumerate(queries_to_run):
                skip = 0
                total_found = 0
                query_results = 0
                
                # Build final query with date filter for daily mode
                if mode == "daily":
                    week_ago = (datetime.utcnow() - timedelta(days=7)).strftime("%Y%m%d")
                    final_query = f"{query}+AND+effective_time:[{week_ago}+TO+*]"
                else:
                    final_query = query

                # 1. Check Total Count first (Smart Check) - with retry
                check_params = {"search": final_query, "limit": 1}
                check_res = await self._fetch_with_retry(client, self.BASE_URL, check_params, self.MAX_RETRIES)
                
                if check_res is None:
                    continue
                    
                if check_res.status_code == 200:
                    meta = check_res.json().get("meta", {})
                    total_found = meta.get("results", {}).get("total", 0)
                    if total_found > 0:
                        logger.info(f"🔍 [{idx+1}/{total_queries}] '{query[:40]}...' → {total_found} records")
                elif check_res.status_code == 404:
                    continue
                else:
                    logger.warning(f"OpenFDA Check failed for '{query[:40]}': {check_res.status_code}")
                    continue

                if total_found == 0:
                    continue

                # Rate limiting - 쿼리 사이에 대기
                await asyncio.sleep(self.REQUEST_DELAY)

                # 2. Fetch Loop (All Pages) - with retry
                max_skip = 1000  # OpenFDA has a max skip limit
                while skip < min(total_found, max_skip):
                    params = {
                        "search": final_query,
                        "limit": min(limit, 100),  # Max 100 per request
                        "skip": skip
                    }
                    
                    res = await self._fetch_with_retry(client, self.BASE_URL, params, self.MAX_RETRIES)
                    
                    if res is None:
                        break
                    
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
                            break  # End of results
                        
                        skip += limit
                        await asyncio.sleep(self.REQUEST_DELAY)  # Rate limiting - 1초 대기
                    elif res.status_code == 404:
                        break
                    else:
                        logger.error(f"OpenFDA Fetch Error: {res.status_code}")
                        break
                
                if query_results > 0:
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
