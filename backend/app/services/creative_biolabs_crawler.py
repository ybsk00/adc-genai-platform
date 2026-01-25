import asyncio
import random
import logging
import json
import time
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

logger = logging.getLogger(__name__)

class CreativeBiolabsCrawler:
    """
    Creative Biolabs ADC/Antibody Stealth Crawler & AI Enrichment (V3 Fixed)
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
        self.global_semaphore = asyncio.Semaphore(2)
        genai.configure(api_key=settings.GOOGLE_API_KEY)
        self.model = genai.GenerativeModel('gemini-1.5-flash')

    async def _init_browser(self, p) -> BrowserContext:
        browser = await p.chromium.launch(headless=True, args=['--no-sandbox'])
        context = await browser.new_context(user_agent=self.ua.random)
        
        async def route_intercept(route):
            if route.request.resource_type in ["image", "media", "font", "stylesheet"]:
                await route.abort()
            else:
                await route.continue_()
        await context.route("**/*", route_intercept)
        return context

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10) -> int:
        logger.info(f"🚀 [Start] {category_name}")
        count = 0
        async with async_playwright() as p:
            context = await self._init_browser(p)
            page = await context.new_page()
            try:
                await page.goto(base_url, wait_until="domcontentloaded", timeout=30000)
                links = await page.evaluate("""
                    () => Array.from(document.querySelectorAll('a[href*="/adc/"], a[href*="/p-"]'))
                        .map(a => a.href)
                        .filter(href => href.includes('.htm') && !href.includes('index.htm') && !href.includes('classify-'))
                """)
                links = list(set(links))[:limit]
                for link in links:
                    res = await self._process_single_product(context, link, category_name)
                    if res:
                        await self._enrich_and_save_single(res)
                        count += 1
            except Exception as e:
                logger.error(f"Error in crawl_category: {e}")
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
                specs = await page.evaluate("""
                    () => {
                        const data = {};
                        const tables = Array.from(document.querySelectorAll('table'));
                        tables.forEach(table => {
                            const rows = Array.from(table.querySelectorAll('tr'));
                            rows.forEach(row => {
                                const cells = Array.from(row.querySelectorAll('td, th'));
                                if (cells.length >= 2) {
                                    const key = cells[0].innerText.replace(':', '').trim();
                                    const val = cells[1].innerText.trim();
                                    if (key.length < 50) data[key] = val;
                                }
                            });
                        });
                        return data;
                    }
                """)
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
            except: return None
            finally: await page.close()

    async def _enrich_with_gemini(self, raw_data: Dict) -> Dict:
        """AI Deep Enrichment - Extracting Biological Indicators"""
        text = f"Product: {raw_data['product_name']}\nSpecs: {json.dumps(raw_data['specs'])}\nText: {raw_data['body_text'][:2000]}"
        
        prompt = """
        Extract ADC metadata from text.
        JSON format:
        {
            "target": "Antigen Symbol",
            "properties": {
                "binding_affinity": "Kd value",
                "isotype": "...",
                "host": "..."
            },
            "summary": "1-sentence summary"
        }
        Text: """ + text
        
        try:
            response = await self.model.generate_content_async(
                prompt,
                generation_config=genai.GenerationConfig(response_mime_type="application/json")
            )
            return json.loads(response.text)
        except:
            return {}

    async def _enrich_and_save_single(self, raw_data):
        try:
            ai_data = await self._enrich_with_gemini(raw_data)
            
            final_data = {
                "ambeed_cat_no": raw_data["ambeed_cat_no"],
                "product_name": raw_data["product_name"],
                "product_url": raw_data["product_url"],
                "category": raw_data["category"],
                "source_name": raw_data["source_name"],
                "target": ai_data.get("target"),
                "summary": ai_data.get("summary"),
                "properties": {
                    **raw_data.get("specs", {}),
                    **ai_data.get("properties", {})
                },
                "crawled_at": raw_data["crawled_at"]
            }
            
            # Embedding
            embed_text = f"{final_data['product_name']} {final_data.get('target', '')} {final_data.get('summary', '')}"
            final_data["embedding"] = await rag_service.generate_embedding(embed_text)

            supabase.table("commercial_reagents").upsert(final_data, on_conflict="ambeed_cat_no").execute()
            logger.info(f"      ✅ Saved: {final_data['product_name']} (Target: {final_data.get('target')})")
        except Exception as e:
            logger.error(f"      ❌ Save failed: {e}")

    async def run(self, search_term: str, limit: int, job_id: str):
        from app.api.scheduler import update_job_status
        await update_job_status(job_id, status="running")
        try:
            targets = {}
            if search_term and search_term != 'all':
                for cat, url in self.CATEGORIES.items():
                    if search_term.lower() in cat.lower():
                        targets[cat] = url
            else:
                targets = self.CATEGORIES

            results = []
            for cat, url in targets.items():
                res = await self.crawl_category(cat, url, limit)
                results.append(res)
                
            total = sum(results)
            await update_job_status(job_id, status="completed", records_drafted=total, completed_at=datetime.utcnow().isoformat())
        except Exception as e:
            logger.error(f"Run failed: {e}")
            await update_job_status(job_id, status="failed", errors=[str(e)])

creative_crawler = CreativeBiolabsCrawler()