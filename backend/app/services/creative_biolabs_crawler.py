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
from app.services.job_lock import job_lock

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
        # genai.configure is NOT called here to avoid fork/gRPC issues
        self.model = None

    def _get_model(self):
        """Lazy initialization of Gemini model to prevent gRPC fork issues"""
        if not self.model:
            genai.configure(api_key=settings.GOOGLE_API_KEY)
            self.model = genai.GenerativeModel('gemini-1.5-flash')
        return self.model

    async def _init_browser(self, p) -> BrowserContext:
        try:
            logger.info("🌐 Launching Browser (Headless)...")
            browser = await p.chromium.launch(
                headless=True, 
                args=[
                    '--no-sandbox', 
                    '--disable-setuid-sandbox', 
                    '--disable-dev-shm-usage',
                    '--disable-gpu',
                    '--no-zygote',
                    '--single-process'
                ],
                timeout=60000
            )
            context = await browser.new_context(user_agent=self.ua.random)
            
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

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10) -> int:
        logger.info(f"🚀 [Start] {category_name}")
        count = 0
        async with async_playwright() as p:
            try:
                context = await self._init_browser(p)
            except Exception:
                return 0

            page = await context.new_page()
            try:
                logger.info(f"📄 Navigating to {base_url}...")
                await page.goto(base_url, wait_until="domcontentloaded", timeout=30000)
                links = await page.evaluate("""
                    () => Array.from(document.querySelectorAll('a[href*="/adc/"], a[href*="/p-"]'))
                        .map(a => a.href)
                        .filter(href => href.includes('.htm') && !href.includes('index.htm') && !href.includes('classify-'))
                """)
                links = list(set(links))[:limit]
                logger.debug(f"🔎 Found {len(links)} links in {category_name}")
                if not links:
                    logger.warning(f"⚠️ No links found for {category_name} at {base_url}")
                for link in links:
                    res = await self._process_single_product(context, link, category_name)
                    if res:
                        await self._enrich_and_save_single(res)
                        count += 1
            except Exception as e:
                logger.error(f"Error in crawl_category {category_name}: {e}", exc_info=True)
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
            except Exception as e:
                logger.error(f"Failed to process {url}: {e}")
                return None
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
            model = self._get_model()
            response = await model.generate_content_async(
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

            res = supabase.table("commercial_reagents").upsert(final_data, on_conflict="ambeed_cat_no").execute()
            
            if res.data:
                logger.info(f"[DB_SAVE_SUCCESS] {final_data['product_name']} (Target: {final_data.get('target')})")
                # [실시간 AI 정제 트리거]
                asyncio.create_task(self._trigger_refinement(res.data[0]))
            
        except Exception as e:
            logger.error(f"      ❌ Save failed: {e}")

    async def _trigger_refinement(self, record: Dict):
        """AI 정제 엔진 비동기 호출 (상세 로그 추가)"""
        try:
            logger.debug(f"🔍 [AI Trigger] Starting refinement for: {record.get('product_name')} (ID: {record.get('id')})")
            
            # AI 분석 서비스 호출
            analysis = await ai_refiner.refine_single_record(record)
            
            if analysis and "error" not in analysis:
                update_data = {
                    "target": analysis.get("target"),
                    "ai_refined": True,
                    "properties": {**record.get("properties", {}), "ai_analysis": analysis}
                }
                res = supabase.table("commercial_reagents").update(update_data).eq("id", record["id"]).execute()
                if res.data:
                    logger.info(f"✨ [AI Success] Refined product: {record.get('product_name')}")
                else:
                    logger.error(f"❌ [AI DB Error] Failed to update record: {record.get('id')}")
            else:
                error_msg = analysis.get("error") if analysis else "Empty response"
                logger.warning(f"⚠️ [AI Skip/Fail] Refinement skipped for {record.get('product_name')}: {error_msg}")
                
        except Exception as e:
            logger.error(f"🔥 [AI Critical Error] Background refinement failed for {record.get('id')}: {str(e)}", exc_info=True)

    async def run(self, search_term: str, limit: int, job_id: str):
        """Orchestrate the crawling process with Stability & Hanging Prevention"""
        from app.api.scheduler import update_job_status
        await update_job_status(job_id, status="running")
        
        logger.info(f"🚀 [CRAWLER START] Job: {job_id} | Term: {search_term} | Limit: {limit}")
        
        # --- Sequential Execution Lock ---
        lock_acquired = False
        wait_start = time.time()
        while not lock_acquired:
            if await job_lock.acquire("global_crawler_lock"):
                lock_acquired = True
            else:
                if time.time() - wait_start > 3600: # 1 hour timeout
                    await update_job_status(job_id, status="failed", errors=["Timed out waiting for global crawler lock"])
                    return
                logger.info(f"⏳ Job {job_id} waiting for global crawler lock...")
                await asyncio.sleep(30)
        
        total_processed = 0
        
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
                try:
                    # Give each category a max of 15 minutes
                    res = await asyncio.wait_for(self.crawl_category(cat, url, limit), timeout=900)
                    results.append(res)
                except Exception as cat_e:
                    logger.error(f"❌ Category {cat} failed or timed out: {cat_e}")
                    results.append(0)
                
            total_processed = sum(results)
            logger.info(f"🏁 [CRAWLER FINISHED] Total Records: {total_processed}")
            
            await update_job_status(
                job_id, 
                status="completed", 
                records_drafted=total_processed, 
                completed_at=datetime.utcnow().isoformat()
            )
        except Exception as e:
            logger.error(f"🔥 Crawler run failed: {e}", exc_info=True)
            await update_job_status(job_id, status="failed", errors=[str(e)])
        finally:
            # Release Lock
            await job_lock.release("global_crawler_lock")
            
            # Safety Fallback
            job = supabase.table("sync_jobs").select("status").eq("id", job_id).execute()
            if job.data and job.data[0]['status'] == 'running':
                await update_job_status(job_id, status="completed", message="Process ended (Fallback)")

creative_crawler = CreativeBiolabsCrawler()