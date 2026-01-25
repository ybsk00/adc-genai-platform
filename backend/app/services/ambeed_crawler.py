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

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10) -> int:
        logger.info(f"🚀 [AMBEED SMILES CRAWL] {category_name}")
        count = 0
        page_num = 1
        async with async_playwright() as p:
            try:
                context = await self._init_browser(p)
                page = await context.new_page()
                while count < limit:
                    separator = "&" if "?" in base_url else "?"
                    url = base_url if page_num == 1 else f"{base_url}{separator}page={page_num}"
                    logger.info(f"📂 Navigating to Page {page_num}...")
                    try:
                        await page.goto(url, wait_until="domcontentloaded", timeout=45000)
                        await asyncio.sleep(2)
                        links = await page.evaluate("""
                            () => Array.from(document.querySelectorAll('a[href*="/products/"], a[href*="/record/"]'))
                                .map(a => a.href)
                                .filter(href => !href.includes('javascript') && !href.includes('google'))
                        """)
                        product_links = list(set(links))
                        if not product_links: break
                        for link in product_links:
                            if count >= limit: break
                            res = await self._process_single_product(context, link, category_name)
                            if res:
                                await self._enrich_and_save_single(res)
                                count += 1
                        page_num += 1
                    except: break
            finally: await context.close()
        return count

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

    async def _enrich_and_save_single(self, raw_data):
        try:
            # 1. SMILES 검증 및 보완
            smiles = raw_data.get("smiles_code")
            is_valid = self.validate_smiles(smiles)
            
            if not is_valid:
                logger.warning(f"⚠️ Invalid or Missing SMILES for {raw_data['product_name']}. Trying PubChem...")
                fallback_smiles = await self.fetch_pubchem_smiles_by_cas(raw_data.get("cas_number"))
                if fallback_smiles:
                    smiles = fallback_smiles
                    is_valid = True
                else:
                    logger.error(f"❌ Failed to obtain valid SMILES for {raw_data['product_name']}")
            
            # 2. AI 정제 (IC50 등)
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
            
            # 3. 임베딩 및 저장
            embed_text = f"{final_data['product_name']} {final_data.get('smiles_code') or ''}"
            final_data["embedding"] = await rag_service.generate_embedding(embed_text)
            supabase.table("commercial_reagents").upsert(final_data, on_conflict="ambeed_cat_no").execute()
            
            status = "✅ [SMILES OK]" if is_valid else "⚠️ [SMILES MISSING]"
            logger.info(f"{status} {final_data['product_name']}")
            
        except Exception as e:
            logger.error(f"Save failed: {e}")

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
            total += await self.crawl_category(cat, url, limit)
        await update_job_status(job_id, status="completed", records_drafted=total, completed_at=datetime.utcnow().isoformat())

ambeed_crawler = AmbeedCrawler()