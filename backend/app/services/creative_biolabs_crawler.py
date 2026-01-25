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

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 1000, job_id: str = None) -> int:
        logger.info(f"🚀 [CRAWL START] Category: {category_name} | URL: {base_url}")
        from app.api.scheduler import update_job_status, get_job_from_db, is_cancelled
        
        count = 0
        batch_data = []
        start_page_offset = 0
        
        # 1. Offset 관리
        if job_id:
            job_data = await get_job_from_db(job_id)
            if job_data and job_data.get("last_processed_page"):
                start_page_offset = job_data["last_processed_page"]
                logger.info(f"⏭️ Skipping {start_page_offset} pages based on last_processed_page")

        async with async_playwright() as p:
            try:
                context = await self._init_browser(p)
                page = await context.new_page()
                current_url = base_url.strip().replace("biollabs", "biolabs").replace("--", "-").replace("..", ".")
                processed_pages_count = 0
                
                while current_url and count < limit:
                    if job_id and await is_cancelled(job_id):
                        logger.info("🛑 Job cancelled.")
                        break

                    # Skip pages based on offset
                    if processed_pages_count < start_page_offset:
                        logger.info(f"⏩ Skipping page (offset): {current_url}")
                        await page.goto(current_url, wait_until="domcontentloaded", timeout=45000)
                        next_page = await self._get_next_page_url(page)
                        current_url = next_page if next_page and next_page != current_url else None
                        processed_pages_count += 1
                        continue

                    logger.info(f"📄 Scraping Page: {current_url} (Page Index: {processed_pages_count})")
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
                            
                            # 2. ID 기반 필터링
                            cat_no = f"CB-{link.split('/')[-1].replace('.htm', '')}"
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
                                    
                                    # 3. 2개 단위 배치 저장 - 사장님 지시 하향 조정
                                    if len(batch_data) >= 2:
                                        await self._save_batch(batch_data)
                                        logger.info(f"📢 [보고] CB 배치 저장 완료. 현재 총 수집: {count}")
                                        batch_data = []
                                        if job_id:
                                            await update_job_status(job_id, records_drafted=count, last_processed_page=processed_pages_count)

                            await asyncio.sleep(0.5)
                        
                        processed_pages_count += 1
                        if job_id:
                            await update_job_status(job_id, last_processed_page=processed_pages_count)

                        next_page = await self._get_next_page_url(page)
                        current_url = next_page if next_page and next_page != current_url else None
                        
                    except Exception as page_e:
                        logger.error(f"Page processing error: {page_e}")
                        break
                
                # 남은 데이터 저장
                if batch_data:
                    await self._save_batch(batch_data)
                    if job_id:
                        await update_job_status(job_id, records_drafted=count)
            except Exception as e:
                logger.error(f"Category crawl error: {e}", exc_info=True)
            finally:
                await context.close()
        return count

    async def _get_next_page_url(self, page):
        return await page.evaluate(r"""
            () => {
                const nextBtn = Array.from(document.querySelectorAll('a')).find(a => 
                    ['Next', '>', 'Next Page'].includes(a.innerText.trim())
                );
                return nextBtn ? nextBtn.href : null;
            }
        """
        )

    async def _enrich_and_prepare_item(self, raw_data):
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
            return final_data
        except Exception as e:
            logger.error(f"Prepare item failed: {e}")
            return None

    async def _save_batch(self, items: List[Dict]):
        try:
            if not items: return
            supabase.table("commercial_reagents").upsert(items, on_conflict="ambeed_cat_no").execute()
        except Exception as e:
            logger.error(f"Batch save failed: {e}")

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

    async def run(self, search_term: str, limit: int, job_id: str):
        from app.api.scheduler import update_job_status
        await update_job_status(job_id, status="running")
        targets = {}
        if search_term and search_term != 'all':
            for cat, url in self.CATEGORIES.items():
                if search_term.lower() in cat.lower(): targets[cat] = url
        else: targets = self.CATEGORIES
        total = 0
        for cat, url in targets.items(): total += await self.crawl_category(cat, url, limit, job_id)
        await update_job_status(job_id, status="completed", records_drafted=total, completed_at=datetime.utcnow().isoformat())

creative_crawler = CreativeBiolabsCrawler()
