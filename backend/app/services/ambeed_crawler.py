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
            self.model = genai.GenerativeModel('gemini-1.5-flash')
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

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10, job_id: str = None) -> int:
        logger.info(f"🚀 [AMBEED SMILES CRAWL] {category_name}")
        from app.api.scheduler import update_job_status, get_job_from_db, is_cancelled
        
        count = 0
        start_page = 1
        
        # 1. Offset 관리: DB에서 마지막 진행 페이지 조회
        if job_id:
            job_data = await get_job_from_db(job_id)
            if job_data and job_data.get("last_processed_page"):
                start_page = job_data["last_processed_page"] + 1
                logger.info(f"⏭️ Resuming from page {start_page}")

        page_num = start_page
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
                    logger.info(f"📂 Navigating to Page {page_num}... (Current Count: {count}/{limit})")
                    
                    try:
                        await page.goto(url, wait_until="domcontentloaded", timeout=45000)
                        await asyncio.sleep(2)
                        
                        # 제품 링크 및 고유 번호 동시 추출 시도 (필터링 효율화)
                        products = await page.evaluate("""
                            () => {
                                return Array.from(document.querySelectorAll('.product-item, .item')).map(el => {
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
                            # Fallback for simple links if structured parsing fails
                            links = await page.evaluate("""
                                () => Array.from(document.querySelectorAll('a[href*="/products/"], a[href*="/record/"]'))
                                    .map(a => a.href)
                                    .filter(href => !href.includes('javascript') && !href.includes('google'))
                            """)
                            products = [{"href": link, "cat_no": None} for link in set(links)]

                        if not products:
                            logger.info(f"🏁 No more products found on page {page_num}. Ending.")
                            break
                        
                        for prod in products:
                            if count >= limit: break
                            
                            link = prod["href"]
                            cat_no = prod["cat_no"] or link.split('/')[-1].replace('.html', '')
                            
                            # 2. ID 기반 필터링: 이미 DB에 존재하면 Skip
                            existing = supabase.table("commercial_reagents").select("id").eq("ambeed_cat_no", cat_no).execute()
                            if existing.data:
                                logger.info(f"⏩ Skipping existing product: {cat_no}")
                                continue

                            res = await self._process_single_product(context, link, category_name)
                            if res:
                                final_item = await self._enrich_and_prepare_item(res)
                                if final_item:
                                    batch_data.append(final_item)
                                    count += 1
                                    
                                    # 3. 20개 단위 배치 저장 (Batch Upsert)
                                    if len(batch_data) >= 20:
                                        await self._save_batch(batch_data)
                                        batch_data = [] # Clear memory
                                        if job_id:
                                            await update_job_status(job_id, records_drafted=count, last_processed_page=page_num)
                                            logger.info(f"💾 20 items batch saved. Current count: {count}")

                        # 페이지 종료 후 상태 업데이트
                        if job_id:
                            await update_job_status(job_id, last_processed_page=page_num)
                        
                        page_num += 1
                    except Exception as e:
                        logger.error(f"❌ Error on page {page_num}: {e}")
                        break
                
                # 남은 데이터 저장
                if batch_data:
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
            
            # 2. AI 정제 (Gemini)
            ai_data = await self._enrich_with_gemini(raw_data)
            
            final_data = {
                "ambeed_cat_no": raw_data["ambeed_cat_no"],
                "product_name": raw_data["product_name"],
                "product_url": raw_data["product_url"],
                "category": raw_data["category"],
                "source_name": "Ambeed",
                "smiles_code": smiles if is_valid else None,
                "cas_number": raw_data.get("cas_number"),
                "target": ai_data.get("target"),
                "summary": ai_data.get("summary"),
                "properties": ai_data.get("properties", {}),
                "crawled_at": raw_data["crawled_at"]
            }
            
            # 임베딩 생성
            embed_text = f"{final_data['product_name']} {final_data.get('smiles_code') or ''} {final_data.get('target') or ''}"
            final_data["embedding"] = await rag_service.generate_embedding(embed_text)
            
            return final_data
        except Exception as e:
            logger.error(f"Failed to prepare item {raw_data.get('ambeed_cat_no')}: {e}")
            return None

    async def _save_batch(self, items: List[Dict]):
        """배치 UPSERT 실행"""
        try:
            if not items: return
            supabase.table("commercial_reagents").upsert(items, on_conflict="ambeed_cat_no").execute()
        except Exception as e:
            logger.error(f"Batch save failed: {e}")

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

    async def _enrich_with_gemini(self, raw_data: Dict) -> Dict:
        prompt = f"Extract target, ic50/solubility, and summary from: {raw_data['body_text'][:1500]}\nJSON Output: {{ 'target': '...', 'properties': {{ 'ic50': '...', 'solubility': '...' }}, 'summary': '...' }}"
        try:
            model = self._get_model()
            response = await model.generate_content_async(prompt, generation_config=genai.GenerationConfig(response_mime_type="application/json"))
            return json.loads(response.text)
        except: return {}

    async def run(self, search_term: str, limit: int, job_id: str):
        from app.api.scheduler import update_job_status
        await update_job_status(job_id, status="running")
        targets = {cat: url for cat, url in self.CATEGORIES.items() if not search_term or search_term == 'all' or search_term.lower() in cat.lower()}
        total = 0
        for cat, url in targets.items():
            total += await self.crawl_category(cat, url, limit, job_id)
        await update_job_status(job_id, status="completed", records_drafted=total, completed_at=datetime.utcnow().isoformat())

ambeed_crawler = AmbeedCrawler()