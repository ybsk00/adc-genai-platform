import asyncio
import random
import logging
import json
import time
import subprocess
import os

# Force gRPC to be less aggressive
os.environ["GRPC_TYPE_CHECK"] = "0"
from typing import Dict, List, Optional, Any
from datetime import datetime

from playwright.async_api import async_playwright, Page, BrowserContext, TimeoutError as PlaywrightTimeoutError
from fake_useragent import UserAgent
import aiohttp
import google.generativeai as genai

from app.core.config import settings
from app.core.supabase import supabase
from app.services.ai_refiner import ai_refiner
from app.services.rag_service import rag_service
from app.services.job_lock import job_lock

logger = logging.getLogger(__name__)

class AmbeedCrawler:
    """
    Ambeed Full-Scale Crawler (V3 Enhanced)
    Features: Stable Pagination, IC50/Solubility/Cleavability extraction via AI
    """
    
    BASE_URL = "https://www.ambeed.com"
    
    CATEGORIES = {
        "ADC Toxins": "https://www.ambeed.com/adc-toxins.html",
        "ADC Linker": "https://www.ambeed.com/adc-linkers.html", 
        "ADC Related": "https://www.ambeed.com/antibody-drug-conjugate-adc-related.html",
        "Payload": "https://www.ambeed.com/payloads.html",
        "Payload-Linker": "https://www.ambeed.com/payload-linker-conjugate.html",
        "PROTAC Linker": "https://www.ambeed.com/protac-linker.html"
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

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10) -> int:
        logger.info(f"🚀 [AMBEED START] Category: {category_name}")
        count = 0
        page_num = 1
        
        async with async_playwright() as p:
            try:
                context = await self._init_browser(p)
                page = await context.new_page()
                
                while count < limit:
                    separator = "&" if "?" in base_url else "?"
                    url = base_url if page_num == 1 else f"{base_url}{separator}page={page_num}"
                    
                    logger.info(f"📂 Navigating to {category_name} Page {page_num}...")
                    try:
                        await page.goto(url, wait_until="domcontentloaded", timeout=45000)
                        await asyncio.sleep(2)
                        
                        # 상품 링크 추출
                        links = await page.evaluate("""
                            () => Array.from(document.querySelectorAll('a[href*="/products/"], a[href*="/record/"]'))
                                .map(a => a.href)
                                .filter(href => !href.includes('javascript') && !href.includes('google'))
                        """)
                        product_links = list(set(links))
                        
                        if not product_links:
                            logger.info(f"🏁 No more products on page {page_num}. Ending.")
                            break

                        for link in product_links:
                            if count >= limit: break
                            res = await self._process_single_product(context, link, category_name)
                            if res:
                                await self._enrich_and_save_single(res)
                                count += 1
                            await asyncio.sleep(0.3)

                        page_num += 1
                    except Exception as page_e:
                        logger.error(f"Error on page {page_num}: {page_e}")
                        break
                        
            finally:
                await context.close()
        return count

    async def _process_single_product(self, context, url, category):
        async with self.global_semaphore:
            page = await context.new_page()
            try:
                await page.goto(url, wait_until="domcontentloaded", timeout=30000)
                title = await page.title()
                body = await page.inner_text("body")
                
                extracted = await page.evaluate("""
                    () => {
                        const getText = (label) => {
                            const tds = Array.from(document.querySelectorAll('td'));
                            for (let i = 0; i < tds.length; i++) {
                                if (tds[i].innerText.includes(label)) return tds[i+1]?.innerText.trim();
                            }
                            return null;
                        };
                        return {
                            cat_no: getText('Catalog No') || getText('Cat. No'),
                            cas_no: getText('CAS No'),
                            formula: getText('Formula'),
                            mw: getText('M.W'),
                            smiles: getText('SMILES'),
                        };
                    }
                """)
                
                return {
                    "ambeed_cat_no": extracted['cat_no'] or url.split('/')[-1].replace('.html', ''),
                    "cas_number": extracted['cas_no'],
                    "product_name": title.split('|')[0].strip(),
                    "product_url": url,
                    "category": category,
                    "smiles_code": extracted['smiles'],
                    "properties": {
                        "formula": extracted['formula'],
                        "mw": extracted['mw']
                    },
                    "body_text": body[:2500],
                    "source_name": "Ambeed",
                    "crawled_at": datetime.utcnow().isoformat()
                }
            except Exception as e:
                logger.error(f"Failed to process {url}: {e}")
                return None
            finally: await page.close()

    async def _enrich_with_gemini(self, raw_data: Dict) -> Dict:
        """AI Deep Enrichment - 정량 지표 최우선 추출"""
        cat = raw_data.get("category", "").lower()
        
        # 특화 지표 추출 가이드
        goal = "Extract basic properties."
        if "toxin" in cat or "payload" in cat:
            goal = "CRITICAL: Look for IC50 or EC50 values (potency/toxicity)."
        elif "linker" in cat:
            goal = "CRITICAL: Look for Solubility (in DMSO/Water) and Cleavability (Cleavable/Non-cleavable type)."

        prompt = f"""
        Extract ADC reagent data.
        Requirement: {goal}
        
        Output JSON:
        {{
            "target": "Target Antigen (if any)",
            "properties": {{
                "ic50": "potency value",
                "solubility": "solubility description",
                "cleavability": "cleavage type",
                "mw": "...",
                "purity": "..."
            }},
            "summary": "1-sentence description"
        }}
        Text: {raw_data['body_text']}
        """
        
        try:
            model = self._get_model()
            response = await model.generate_content_async(prompt, generation_config=genai.GenerationConfig(response_mime_type="application/json"))
            return json.loads(response.text)
        except: return {}

    async def _enrich_and_save_single(self, raw_data):
        try:
            ai_data = await self._enrich_with_gemini(raw_data)
            final_data = {
                "ambeed_cat_no": raw_data["ambeed_cat_no"],
                "product_name": raw_data["product_name"],
                "product_url": raw_data["product_url"],
                "category": raw_data["category"],
                "source_name": "Ambeed",
                "target": ai_data.get("target"),
                "summary": ai_data.get("summary"),
                "properties": {**raw_data.get("properties", {}), **ai_data.get("properties", {})},
                "smiles_code": raw_data.get("smiles_code"),
                "cas_number": raw_data.get("cas_number"),
                "crawled_at": raw_data["crawled_at"]
            }
            embed_text = f"{final_data['product_name']} {final_data.get('target') or ''}"
            final_data["embedding"] = await rag_service.generate_embedding(embed_text)
            supabase.table("commercial_reagents").upsert(final_data, on_conflict="ambeed_cat_no").execute()
            logger.info(f"✅ [AMBEED SAVE] {final_data['product_name']}")
        except Exception as e:
            logger.error(f"Save failed: {e}")

    async def run(self, search_term: str, limit: int, job_id: str):
        from app.api.scheduler import update_job_status
        await update_job_status(job_id, status="running")
        
        targets = {}
        if search_term and search_term != 'all':
            for cat, url in self.CATEGORIES.items():
                if search_term.lower() in cat.lower(): targets[cat] = url
        else: targets = self.CATEGORIES

        total = 0
        for cat, url in targets.items():
            res = await self.crawl_category(cat, url, limit)
            total += res
        
        await update_job_status(job_id, status="completed", records_drafted=total, completed_at=datetime.utcnow().isoformat())

ambeed_crawler = AmbeedCrawler()
