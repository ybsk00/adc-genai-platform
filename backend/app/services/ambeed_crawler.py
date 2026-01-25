import asyncio
import random
import logging
import json
import time
import subprocess
import os
import re

# Force gRPC to be less aggressive
os.environ["GRPC_TYPE_CHECK"] = "0"
from typing import Dict, List, Optional, Any
from datetime import datetime

from playwright.async_api import async_playwright, Page, BrowserContext, TimeoutError as PlaywrightTimeoutError
from fake_useragent import UserAgent
import aiohttp
import google.generativeai as genai
from rdkit import Chem # SMILES Validation

from app.core.config import settings
from app.core.supabase import supabase
from app.services.ai_refiner import ai_refiner
from app.services.rag_service import rag_service
from app.services.job_lock import job_lock

logger = logging.getLogger(__name__)

class AmbeedCrawler:
    """
    Ambeed Full-Scale Crawler (V4 SMILES Specialized)
    - Precise Detail Page Parsing for SMILES
    - RDKit Validation
    - PubChem Fallback (via CAS No)
    """
    
    BASE_URL = "https://www.ambeed.com"
    CATEGORIES = {
        "ADC Toxins": "https://www.ambeed.com/adc-toxins.html",
        "ADC Linker": "https://www.ambeed.com/adc-linkers.html"
    }
    
    def __init__(self):
        self.ua = UserAgent()
        self.global_semaphore = asyncio.Semaphore(1)
        self.model = None

    def _get_model(self):
        if not self.model:
            genai.configure(api_key=settings.GOOGLE_API_KEY)
            # 사장님 지시: 최신 2.0 Flash 모델로 업그레이드
            self.model = genai.GenerativeModel('gemini-2.0-flash')
        return self.model

    async def _init_browser(self, p) -> BrowserContext:
        try:
            browser = await p.chromium.launch(
                headless=True,
                args=['--no-sandbox', '--disable-setuid-sandbox', '--disable-dev-shm-usage', '--disable-gpu']
            )
            context = await browser.new_context(user_agent=self.ua.random)
            
            # --- RESOURCE BLOCKING (Critical for Speed) ---
            async def route_intercept(route):
                if route.request.resource_type in ["image", "media", "font", "stylesheet"]:
                    await route.abort()
                else:
                    await route.continue_()

            await context.route("**/*", route_intercept)
            return context
        except Exception as e:
            logger.error(f"🔥 Ambeed Browser Launch Failed: {e}")
            raise e

    def validate_smiles(self, smiles: str) -> bool:
        """RDKit을 사용한 SMILES 유효성 검증"""
        if not smiles or not isinstance(smiles, str): return False
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol is not None
        except: 
            return False

    async def fetch_pubchem_smiles_by_cas(self, cas: str) -> Optional[str]:
        """CAS No를 이용해 PubChem에서 SMILES 추론 (Fallback)"""
        if not cas: return None
        logger.info(f"🧬 Attempting PubChem fallback for CAS: {cas}")
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas}/property/IsomericSMILES/JSON"
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=10) as response:
                    if response.status == 200:
                        data = await response.json()
                        smiles = data['PropertyTable']['Properties'][0].get('IsomericSMILES')
                        if self.validate_smiles(smiles):
                            logger.info(f"✨ Successfully found SMILES on PubChem for {cas}")
                            return smiles
        except Exception as e:
            logger.warning(f"⚠️ PubChem fallback failed for {cas}: {e}")
        return None

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10, job_id: str = None, start_page: int = 1, batch_size: int = 2) -> int:
        logger.info(f"🚀 [AMBEED SMILES CRAWL] {category_name} (Start Page: {start_page}, Batch: {batch_size})")
        from app.api.scheduler import update_job_status, get_job_from_db, is_cancelled
        
        count = 0
        current_start_page = start_page
        
        # 1. Offset 관리: DB에서 마지막 진행 페이지 조회 (start_page가 1일 때만 DB 조회)
        if job_id and current_start_page == 1:
            job_data = await get_job_from_db(job_id)
            if job_data and job_data.get("last_processed_page"):
                current_start_page = job_data["last_processed_page"] + 1
                logger.info(f"⏭️ Resuming from page {current_start_page}")

        page_num = current_start_page
        batch_data = []
        
        async with async_playwright() as p:
            try:
                context = await self._init_browser(p)
                page = await context.new_page()
                
                while count < limit:
                    if job_id and await is_cancelled(job_id):
                        logger.info(f"🛑 Job {job_id} cancelled.")
                        break

                    separator = "&" if "?" in base_url else "?"
                    url = base_url if page_num == 1 else f"{base_url}{separator}page={page_num}"
                    logger.info(f"🌐 [PAGE {page_num}] 다음 페이지 접속 중: {url} (현재 수집량: {count}/{limit})")
                    
                    try:
                        # 페이지 접속 시도 (타임아웃 강화 및 재시도 로직은 생략하되 확실히 대기)
                        response = await page.goto(url, wait_until="domcontentloaded", timeout=60000)
                        await asyncio.sleep(3) # 로딩 대기 시간 충분히 부여
                        
                        # 상품 목록 추출 (다양한 선택자 시도)
                        products = await page.evaluate("""
                            () => {
                                // 1. 표준 제품 아이템 클래스
                                let items = Array.from(document.querySelectorAll('.product-item, .item, .product-info-main, .product-item-info'));
                                
                                // 2. 만약 없다면 모든 제품 상세 링크 찾기
                                if (items.length === 0) {
                                    return Array.from(document.querySelectorAll('a[href*="/products/"]'))
                                        .map(a => ({ href: a.href, cat_no: null }))
                                        .filter((v, i, a) => a.findIndex(t => t.href === v.href) === i); // 중복 제거
                                }

                                return items.map(el => {
                                    const linkEl = el.querySelector('a[href*="/products/"], a[href*="/record/"]');
                                    const catNoEl = el.innerText.match(/Cat No:?\s*([A-Z0-9-]+)/i);
                                    return {
                                        href: linkEl ? linkEl.href : null,
                                        cat_no: catNoEl ? catNoEl[1] : null
                                    };
                                }).filter(p => p.href);
                            }
                        """)
                        
                        if not products:
                            # <a> 태그만이라도 뒤져서 찾기
                            links = await page.evaluate("""
                                () => Array.from(document.querySelectorAll('a[href*="/products/"]'))
                                    .map(a => a.href)
                                    .filter(href => href.includes('.html'))
                            """)
                            products = [{"href": link, "cat_no": None} for link in set(links)]

                        if not products:
                            logger.warning(f"⚠️ [PAGE {page_num}] 상품을 찾지 못했습니다. 1페이지를 더 건너뛰어 봅니다.")
                            page_num += 1
                            if page_num > start_page + 50: # 너무 많이 비어있으면 종료
                                logger.error("🏁 연속된 빈 페이지 발생으로 종료합니다.")
                                break
                            continue
                        
                        logger.info(f"📦 [PAGE {page_num}] {len(products)}개의 상품을 찾았습니다. 수집을 시작합니다.")

                        for prod in products:
                            if count >= limit: break
                            
                            link = prod["href"]
                            cat_no = prod["cat_no"] or link.split('/')[-1].replace('.html', '')
                            
                            logger.info(f"🔄 [PROCESS] Item: {cat_no} (Page {page_num})")

                            res = await self._process_single_product(context, link, category_name)
                            if res:
                                final_item = await self._enrich_and_prepare_item(res)
                                if final_item:
                                    batch_data.append(final_item)
                                    count += 1
                                    
                                    if len(batch_data) >= batch_size:
                                        logger.info(f"💾 [단계 3] {batch_size}개 도달! DB 저장 시도 (Page {page_num})")
                                        await self._save_batch(batch_data)
                                        batch_data = []
                                        if job_id:
                                            await update_job_status(job_id, records_drafted=count, last_processed_page=page_num)
                                else:
                                    logger.warning(f"⚠️ [SKIP] Enrichment failed for {cat_no}")
                            else:
                                logger.warning(f"⚠️ [SKIP] Processing failed for {cat_no}")

                        # 한 페이지 처리가 끝나면 "무조건" 페이지 번호 증가
                        logger.info(f"✅ [PAGE {page_num}] 처리 완료. 다음 페이지({page_num + 1})로 이동합니다.")
                        page_num += 1
                        
                        if job_id:
                            await update_job_status(job_id, last_processed_page=page_num)

                    except Exception as e:
                        logger.error(f"❌ Error on page {page_num}: {e}")
                        page_num += 1 # 에러 나도 다음 페이지 시도
                        await asyncio.sleep(5)
                
                # 남은 데이터 저장
                if batch_data:
                    logger.info(f"💾 [Final Batch] 남은 {len(batch_data)}개 저장 시도")
                    await self._save_batch(batch_data)
                    if job_id:
                        await update_job_status(job_id, records_drafted=count)
            finally: 
                await context.close()
        return count

    async def _enrich_and_prepare_item(self, raw_data):
        """저장 전 데이터 보강 및 검증 (Upsert용 데이터 생성)"""
        try:
            # 1. SMILES 검증 및 보완
            smiles = raw_data.get("smiles_code")
            is_valid = self.validate_smiles(smiles)
            
            if not is_valid:
                logger.warning(f"⚠️ Invalid SMILES for {raw_data['ambeed_cat_no']}. Trying PubChem...")
                fallback_smiles = await self.fetch_pubchem_smiles_by_cas(raw_data.get("cas_number"))
                if fallback_smiles:
                    smiles = fallback_smiles
                    is_valid = True
            
            # SMILES 필수 체크 로그
            if not is_valid:
                logger.error(f"❌ SMILES MISSING for {raw_data['ambeed_cat_no']} after fallback.")
            
            # 2. AI 정제 (Gemini 2.0 Flash) - 3단 SMILES 분리 로직 탑재
            ai_data = await self._enrich_with_gemini(raw_data, smiles)
            
            # DB 컬럼 누락 대비: properties 안에 상세 SMILES 정보 백업
            extended_properties = raw_data.get("properties", {})
            if isinstance(extended_properties, dict):
                extended_properties.update({
                    "payload_smiles": ai_data.get("payload_smiles"),
                    "linker_smiles": ai_data.get("linker_smiles"),
                    "full_smiles": ai_data.get("full_smiles") or smiles
                })

            final_data = {
                "ambeed_cat_no": raw_data["ambeed_cat_no"],
                "product_name": raw_data["product_name"],
                "product_url": raw_data["product_url"],
                "category": raw_data["category"],
                "source_name": "Ambeed",
                "smiles_code": smiles if is_valid else None,
                # 아래 컬럼들은 DB에 추가되어야 작동합니다.
                "payload_smiles": ai_data.get("payload_smiles"),
                "linker_smiles": ai_data.get("linker_smiles"),
                "full_smiles": ai_data.get("full_smiles") or smiles,
                "cas_number": raw_data.get("cas_number"),
                "target": ai_data.get("target"),
                "summary": ai_data.get("summary"),
                "properties": extended_properties, # 백업용으로 JSON 안에도 저장
                "crawled_at": raw_data["crawled_at"],
                "ai_refined": True
            }
            
            # 임베딩 생성 (실패해도 저장은 되어야 함)
            try:
                embed_text = f"{final_data['product_name']} {final_data.get('smiles_code') or ''} {final_data.get('target') or ''}"
                # 768 차원인지 확인 등 RAG 서비스 내부 로직에 의존
                embedding = await rag_service.generate_embedding(embed_text)
                if embedding:
                    final_data["embedding"] = embedding
            except Exception as e:
                logger.warning(f"⚠️ Embedding failed for {final_data['ambeed_cat_no']}, proceeding without it: {e}")
            
            return final_data
        except Exception as e:
            logger.error(f"🔥 [치명적 에러] 데이터 준비 단계 실패 ({raw_data.get('ambeed_cat_no')}): {e}", exc_info=True)
            return None

    async def _save_batch(self, items: List[Dict]):
        """배치 UPSERT 실행"""
        try:
            if not items: return False
            
            # DB 연결 상태 체크 (MockClient 방지)
            if hasattr(supabase, "is_mock") and supabase.is_mock:
                logger.error("🔥 [CRITICAL] Supabase is running as a MOCK client. Check .env file!")
                return False

            logger.info(f"📤 Supabase UPSERT 요청 중... ({len(items)}건)")
            res = supabase.table("commercial_reagents").upsert(items, on_conflict="ambeed_cat_no").execute()
            if res.data:
                logger.info(f"✅ DB 저장 완료! ({len(res.data)}건 반영됨)")
                return True
            else:
                logger.error("❌ DB 저장 실패: 응답 데이터가 없습니다.")
                return False
        except Exception as e:
            logger.error(f"🔥 [DB 치명적 에러] Batch save failed: {e}", exc_info=True)
            return False

    async def _enrich_and_save_single(self, raw_data):
        # This is now handled by _enrich_and_prepare_item and _save_batch
        pass

    async def _process_single_product(self, context, url, category):
        async with self.global_semaphore:
            page = await context.new_page()
            try:
                await page.goto(url, wait_until="domcontentloaded", timeout=30000)
                title = await page.title()
                
                # 상세 페이지 정밀 파싱 (SMILES 라벨 타겟팅)
                extracted = await page.evaluate("""
                    () => {
                        const results = { cat_no: '', cas_no: '', smiles: '' };
                        const tds = Array.from(document.querySelectorAll('td, th, span, div'));
                        
                        // SMILES 찾기
                        for (let i = 0; i < tds.length; i++) {
                            const text = tds[i].innerText;
                            if (text.includes('SMILES')) {
                                // 다음 형제 요소나 부모의 텍스트에서 SMILES 추출 시도
                                let smilesVal = tds[i+1]?.innerText || tds[i].parentElement.innerText.split('SMILES')[1];
                                if (smilesVal) results.smiles = smilesVal.replace(':', '').trim();
                            }
                            if (text.includes('CAS No')) {
                                let casVal = tds[i+1]?.innerText || tds[i].parentElement.innerText.split('CAS No')[1];
                                if (casVal) results.cas_no = casVal.replace(':', '').trim();
                            }
                            if (text.includes('Catalog No')) {
                                let catVal = tds[i+1]?.innerText || tds[i].parentElement.innerText.split('Catalog No')[1];
                                if (catVal) results.cat_no = catVal.replace(':', '').trim();
                            }
                        }
                        return results;
                    }
                """)
                
                # 정규식으로 SMILES 한 번 더 정제 (영문/숫자/=/#/()/[]/+/- 등 특수문자 조합)
                raw_smiles = extracted.get('smiles', '')
                if raw_smiles:
                    match = re.search(r'([A-Za-z0-9#%()=+./@\[\]\\-]+)', raw_smiles)
                    if match: raw_smiles = match.group(1)

                return {
                    "ambeed_cat_no": extracted['cat_no'] or url.split('/')[-1].replace('.html', ''),
                    "cas_number": extracted['cas_no'],
                    "product_name": title.split('|')[0].strip(),
                    "product_url": url,
                    "category": category,
                    "smiles_code": raw_smiles,
                    "body_text": await page.inner_text("body"),
                    "source_name": "Ambeed",
                    "crawled_at": datetime.utcnow().isoformat()
                }
            except Exception as e:
                logger.error(f"Failed to process {url}: {e}")
                return None
            finally: await page.close()

    async def _enrich_with_gemini(self, raw_data: Dict, smiles: Optional[str] = None) -> Dict:
        """Gemini를 이용한 3단 SMILES 분리 및 데이터 정제"""
        prompt = f"""
        Extract detailed ADC reagent information from the following content:
        Title: {raw_data['product_name']}
        Body: {raw_data['body_text'][:2000]}
        Input SMILES: {smiles or 'None'}

        Task:
        1. Identify the 'payload_smiles' (the cytotoxic drug part).
        2. Identify the 'linker_smiles' (the chemical linker part).
        3. Identify the 'full_smiles' (the entire drug-linker structure, excluding the antibody).
        4. Extract 'target' and 'summary'.

        Output JSON format:
        {{
            "payload_smiles": "...",
            "linker_smiles": "...",
            "full_smiles": "...",
            "target": "...",
            "summary": "...",
            "properties": {{ "ic50": "...", "solubility": "..." }}
        }}
        If a part is missing, return null. Ensure SMILES strings are valid and complete.
        """
        try:
            model = self._get_model()
            response = await model.generate_content_async(prompt, generation_config=genai.GenerationConfig(response_mime_type="application/json"))
            return json.loads(response.text)
        except Exception as e:
            logger.error(f"Gemini API Error: {e}")
            return {}

    async def run(self, search_term: str, limit: int, job_id: str, start_page: int = 1, batch_size: int = 2):
        from app.api.scheduler import update_job_status
        await update_job_status(job_id, status="running")
        targets = {cat: url for cat, url in self.CATEGORIES.items() if not search_term or search_term == 'all' or search_term.lower() in cat.lower()}
        total = 0
        for cat, url in targets.items():
            total += await self.crawl_category(cat, url, limit, job_id, start_page, batch_size)
        await update_job_status(job_id, status="completed", records_drafted=total, completed_at=datetime.utcnow().isoformat())

ambeed_crawler = AmbeedCrawler()