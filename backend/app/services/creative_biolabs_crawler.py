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

class CreativeBiolabsCrawler:
    """
    Creative Biolabs Full-Scale Crawler (V4 Final Fixed)
    """
    
    BASE_URL = "https://www.creative-biolabs.com/adc/"
    
    CATEGORIES = {
        "ADC Antibody": "https://www.creative-biolabs.com/adc/classify-adc-antibody-products-5.htm",
        "ADC Linker": "https://www.creative-biolabs.com/adc/classify-adc-linker-products-7.htm",
        "ADC Toxin": "https://www.creative-biolabs.com/adc/classify-adc-toxin-products-6.htm",
        "Drug-Linker Complex": "https://www.creative-biolabs.com/adc/classify-drug-linker-complex-products-8.htm"
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
            logger.error(f"🔥 Browser Launch Failed: {e}")
            raise e

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 1000) -> int:
        logger.info(f"🚀 [CRAWL START] Category: {category_name} | URL: {base_url}")
        count = 0
        async with async_playwright() as p:
            try:
                context = await self._init_browser(p)
                page = await context.new_page()
                current_url = base_url.strip().replace("biollabs", "biolabs").replace("--", "-").replace("..", ".")
                processed_pages = set()
                while current_url and count < limit:
                    if current_url in processed_pages: break
                    processed_pages.add(current_url)
                    logger.info(f"📄 Scraping Page: {current_url}")
                    try:
                        await page.goto(current_url, wait_until="domcontentloaded", timeout=45000)
                        await asyncio.sleep(4)
                        links = await page.evaluate(r"""
                            () => Array.from(document.querySelectorAll('a'))
                                .map(a => a.href)
                                .filter(href => {
                                    const isProduct = href.includes('/p-') || (href.includes('/adc/') && href.match(/-[0-9]+\.htm/));
                                    const isClassify = href.includes('classify-');
                                    return isProduct && !isClassify && href.endsWith('.htm');
                                })
                        """
                        )
                        product_links = list(set(links))
                        logger.info(f"🔎 Found {len(product_links)} product links.")
                        for link in product_links:
                            if count >= limit: break
                            res = await self._process_single_product(context, link, category_name)
                            if res:
                                await self._enrich_and_save_single(res)
                                count += 1
                            await asyncio.sleep(0.5)
                        next_page = await page.evaluate(r"""
                            () => {
                                const nextBtn = Array.from(document.querySelectorAll('a')).find(a => 
                                    ['Next', '>', 'Next Page'].includes(a.innerText.trim())
                                );
                                return nextBtn ? nextBtn.href : null;
                            }
                        """
                        )
                        current_url = next_page if next_page and next_page != current_url else None
                    except Exception as page_e:
                        logger.error(f"Page processing error: {page_e}")
                        break
            except Exception as e:
                logger.error(f"Category crawl error: {e}", exc_info=True)
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
                specs = await page.evaluate(r"""
                    () => {
                        const data = {};
                        document.querySelectorAll('tr').forEach(row => {
                            const cells = row.querySelectorAll('td, th');
                            if (cells.length >= 2) {
                                const key = cells[0].innerText.replace(':', '').trim();
                                const val = cells[1].innerText.trim();
                                if (key && val && key.length < 50) data[key] = val;
                            }
                        });
                        return data;
                    }
                """
                )
                return {
                    "ambeed_cat_no": f"CB-{url.split('/')[-1].replace('.htm', '')}",
                    "product_name": title.split('|')[0].strip(),
                    "product_url": url,
                    "category": category,
                    "body_text": body[:3000],
                    "specs": specs,
                    "source_name": "Creative Biolabs",
                    "crawled_at": datetime.utcnow().isoformat()
                }
            except Exception as e:
                logger.error(f"Failed to process {url}: {e}")
                return None
            finally: await page.close()

    async def _enrich_with_gemini(self, raw_data: Dict) -> Dict:
        cat = raw_data.get("category", "").lower()
        goal = "Extract IC50 for Toxins, Solubility/Cleavability for Linkers."
        prompt = f"Extract ADC reagent data. Goal: {goal}\nJSON Output: {{ 'target': '...', 'properties': {{ 'ic50': '...', 'solubility': '...', 'cleavability': '...', 'mw': '...', 'purity': '...' }}, 'summary': '...' }}\nText: {raw_data['body_text'][:2000]}"
        try:
            model = self._get_model()
            response = await model.generate_content_async(prompt, generation_config=genai.GenerationConfig(response_mime_type="application/json"))
            return json.loads(response.text)
        except: return {}

    async def _enrich_and_save_single(self, raw_data):
        try:
            ai_data = await self._enrich_with_gemini(raw_data)
            orig_specs = raw_data.get("specs", {})
            ai_props = ai_data.get("properties", {})
            merged_props = {**orig_specs}
            if isinstance(ai_props, dict): merged_props.update(ai_props)
            final_data = {
                "ambeed_cat_no": raw_data["ambeed_cat_no"],
                "product_name": raw_data["product_name"],
                "product_url": raw_data["product_url"],
                "category": raw_data["category"],
                "source_name": "Creative Biolabs",
                "target": ai_data.get("target"),
                "summary": ai_data.get("summary"),
                "properties": merged_props,
                "crawled_at": raw_data["crawled_at"]
            }
            embed_text = f"{final_data['product_name']} {final_data.get('target') or ''}"
            final_data["embedding"] = await rag_service.generate_embedding(embed_text)
            supabase.table("commercial_reagents").upsert(final_data, on_conflict="ambeed_cat_no").execute()
            logger.info(f"✅ [CB SAVE] {final_data['product_name']}")
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
        for cat, url in targets.items(): total += await self.crawl_category(cat, url, limit)
        await update_job_status(job_id, status="completed", records_drafted=total, completed_at=datetime.utcnow().isoformat())

creative_crawler = CreativeBiolabsCrawler()