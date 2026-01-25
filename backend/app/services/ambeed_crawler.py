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

class AmbeedCrawler:
    """
    Ambeed ADC Reagents Stealth Crawler & AI Enrichment (V2 Enhanced)
    Features: Robust Pagination, Dynamic Scrolling, Concurrency, Anti-blocking
    """
    
    BASE_URL = "https://www.ambeed.com"
    
    # Expanded Categories
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
        self.request_count = 0
        self.break_threshold = 100
        self.batch_size = 5
        self.max_concurrent_categories = 3  # Limit concurrent category crawls
        self.max_concurrent_tabs = 5 # Limit tabs per category or global? Global semaphore is better.
        self.global_semaphore = asyncio.Semaphore(5) # Max 5 concurrent tabs globally

        # Configure Gemini
        genai.configure(api_key=settings.GOOGLE_API_KEY)
        self.model_id = 'gemini-2.5-flash'
        try:
             self.model = genai.GenerativeModel(self.model_id)
        except Exception:
             logger.warning(f"⚠️ Model {self.model_id} not found, falling back to gemini-1.5-flash")
             self.model = genai.GenerativeModel('gemini-1.5-flash')

    async def _init_browser(self, p) -> BrowserContext:
        """Initialize browser with stealth and proxy settings"""
        user_agent = self.ua.random
        
        proxy_config = None
        if settings.PROXY_ENABLED:
            ports = settings.PROXY_PORTS.split(',')
            port = random.choice(ports).strip()
            proxy_config = {
                "server": f"http://{settings.PROXY_HOST}:{port}",
                "username": settings.PROXY_USERNAME,
                "password": settings.PROXY_PASSWORD
            }
            logger.info(f"🌐 Using Proxy: {settings.PROXY_HOST}:{port}")

        # Launch options
        browser = await p.chromium.launch(
            headless=True, # Set False for debugging
            proxy=proxy_config if settings.PROXY_ENABLED else None,
            args=[
                '--disable-blink-features=AutomationControlled',
                '--no-sandbox',
                '--disable-setuid-sandbox',
                '--disable-infobars',
                '--window-position=0,0',
                '--ignore-certificate-errors',
                '--ignore-ssl-errors',
                '--disable-dev-shm-usage'
            ]
        )
        
        context = await browser.new_context(
            user_agent=user_agent,
            viewport={'width': 1920, 'height': 1080},
            locale='en-US',
            timezone_id='America/New_York',
            java_script_enabled=True
        )
        
        # Inject stealth scripts
        await context.add_init_script("""
            Object.defineProperty(navigator, 'webdriver', { get: () => undefined });
            window.navigator.chrome = { runtime: {} };
            Object.defineProperty(navigator, 'languages', { get: () => ['en-US', 'en'] });
            Object.defineProperty(navigator, 'plugins', { get: () => [1, 2, 3, 4, 5] });
        """)
        
        return context

    async def _smart_delay(self, min_s=2.0, max_s=5.0):
        """Random delay to mimic human behavior"""
        delay = random.uniform(min_s, max_s)
        await asyncio.sleep(delay)

    async def _scroll_to_bottom(self, page: Page):
        """Scroll to bottom to trigger lazy loading"""
        try:
            previous_height = await page.evaluate("document.body.scrollHeight")
            while True:
                await page.evaluate("window.scrollTo(0, document.body.scrollHeight)")
                await asyncio.sleep(1.5) # Wait for load
                new_height = await page.evaluate("document.body.scrollHeight")
                if new_height == previous_height:
                    break
                previous_height = new_height
        except Exception as e:
            logger.warning(f"Scroll failed: {e}")

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10) -> int:
        """Crawl a specific category with robust pagination"""
        logger.info(f"🚀 [Start] Category: {category_name}")
        
        count = 0
        page_num = 1
        max_pages = 200 # Safety limit
        is_full_crawl = limit > 1000
        consecutive_empty_pages = 0
        consecutive_load_failures = 0
        
        async with async_playwright() as p:
            context = await self._init_browser(p)
            
            # Create a main listing page
            page = await context.new_page()
            
            try:
                while True:
                    if not is_full_crawl and count >= limit:
                        logger.info(f"🛑 Limit reached for {category_name}. Stopping.")
                        break
                    if page_num > max_pages:
                        logger.info(f"🛑 Max pages reached for {category_name}. Stopping.")
                        break
                    if consecutive_empty_pages >= 3:
                        logger.info(f"🛑 No products found for 3 consecutive pages in {category_name}. Stopping.")
                        break
                    if consecutive_load_failures >= 3:
                        logger.error(f"🛑 Failed to load 3 consecutive pages in {category_name}. Aborting to save resources.")
                        break

                    # Construct URL (Handle ?page= vs &page=)
                    separator = "&" if "?" in base_url else "?"
                    url = base_url if page_num == 1 else f"{base_url}{separator}page={page_num}"
                    
                    logger.info(f"📂 [{category_name}] Navigating to Page {page_num}...")
                    
                    retry_count = 0
                    success = False
                    while retry_count < 3:
                        try:
                            await page.goto(url, wait_until="domcontentloaded", timeout=60000)
                            await self._smart_delay(1, 2)
                            
                            # Check if page is valid (Ambeed sometimes redirects to home on invalid page)
                            if page.url == "https://www.ambeed.com/" and page_num > 1:
                                logger.warning(f"   Redirected to home. Assuming end of pagination.")
                                consecutive_empty_pages = 10 # Force break
                                success = True
                                break

                            # Dynamic Scroll
                            await self._scroll_to_bottom(page)
                            
                            # Wait for product list
                            try:
                                await page.wait_for_selector('.product-item, .pro-list, table.product-table', state='attached', timeout=5000)
                            except:
                                logger.warning(f"   No product list selector found on page {page_num}.")
                            
                            success = True
                            break
                        except Exception as e:
                            logger.warning(f"   Retry {retry_count+1}/3 failed for page {page_num}: {e}")
                            retry_count += 1
                            await asyncio.sleep(5)
                    
                    if not success:
                        logger.error(f"❌ Failed to load page {page_num} after retries. Skipping.")
                        consecutive_load_failures += 1
                        page_num += 1
                        continue
                    
                    # Reset load failure count on success
                    consecutive_load_failures = 0

                    # Extract Links
                    links = await page.evaluate("""
                        Array.from(document.querySelectorAll('a[href*="/products/"], a[href*="/record/"], .product-item a, .pro-list h3 a'))
                            .map(a => a.href)
                            .filter(href => (href.includes('/products/') || href.includes('/record/')) && !href.includes('google') && !href.includes('javascript'))
                    """)
                    links = list(set(links))
                    
                    if not links:
                        logger.info(f"   ⚠️ No links found on page {page_num} ({category_name}).")
                        consecutive_empty_pages += 1
                        page_num += 1
                        continue
                    else:
                        consecutive_empty_pages = 0
                        logger.info(f"   🔎 Found {len(links)} items on page {page_num} ({category_name}). Processing...")

                    # Process Links in Sub-Batches
                    # We use a semaphore to limit concurrent detailed page processing
                    chunk_size = 5
                    for i in range(0, len(links), chunk_size):
                        if not is_full_crawl and count >= limit: break
                        
                        chunk = links[i:i+chunk_size]
                        tasks = []
                        for link in chunk:
                            if not is_full_crawl and count >= limit: break
                            tasks.append(self._process_single_product(context, link, category_name))
                        
                        results = await asyncio.gather(*tasks)
                        valid_results = [r for r in results if r]
                        count += len(valid_results)
                        
                        if valid_results:
                            # Batch Enrich & Save
                            await self._process_batch_enrichment(valid_results)

                    page_num += 1
                    
            except Exception as e:
                logger.error(f"🔥 Critical error in category {category_name}: {e}")
            finally:
                await context.close()
        
        logger.info(f"✅ Finished {category_name}. Total Collected: {count}")
        return count

    async def _process_single_product(self, context: BrowserContext, url: str, category: str) -> Optional[Dict]:
        """Open a new page/tab, extract data, close page."""
        async with self.global_semaphore: # Limit global concurrency
            page = await context.new_page()
            try:
                # logger.debug(f"   Processing: {url}")
                await page.goto(url, wait_until="domcontentloaded", timeout=45000)
                
                # Basic Wait
                try:
                    await page.wait_for_selector('h1', timeout=5000) 
                except:
                    pass

                # Extract Data (Reusing logic but optimized)
                data = await self._extract_html_data_on_page(page, url, category)
                return data

            except Exception as e:
                logger.warning(f"   Failed to process {url}: {e}")
                return None
            finally:
                await page.close()

    async def _extract_html_data_on_page(self, page: Page, url: str, category: str) -> Optional[Dict]:
        """Execute extraction logic on an already open page"""
        try:
            title = await page.title()
            body_text = await page.inner_text("body")
            
            # Efficient extraction using JS evaluation
            extracted = await page.evaluate("""
                () => {
                    const getText = (label) => {
                        // Strategy 1: Table row with label
                        const tds = Array.from(document.querySelectorAll('td'));
                        for (let i = 0; i < tds.length; i++) {
                            if (tds[i].innerText.includes(label)) {
                                if (tds[i+1]) return tds[i+1].innerText.trim();
                            }
                        }
                        // Strategy 2: Strong/B tag followed by text
                        const strongs = Array.from(document.querySelectorAll('strong, b, th'));
                        for (let s of strongs) {
                            if (s.innerText.includes(label)) {
                                let parent = s.parentElement;
                                let text = parent.innerText.replace(s.innerText, '').trim();
                                if (text.startsWith(':')) text = text.substring(1).trim();
                                if (text) return text;
                                // Try next sibling
                                if (s.nextSibling && s.nextSibling.nodeType === 3) return s.nextSibling.nodeValue.trim();
                            }
                        }
                        return null;
                    };

                    return {
                        cat_no: getText('Catalog No') || getText('Cat. No'),
                        cas_no: getText('CAS No') || getText('CAS #'),
                        formula: getText('Molecular Formula') || getText('Formula'),
                        mw: getText('Molecular Weight') || getText('M.W'),
                        smiles: getText('SMILES') || getText('SMILES Code'),
                        purity: getText('Purity'),
                    };
                }
            """)
            
            # Fallbacks
            if not extracted['cat_no']:
                 parts = url.split('/')
                 if parts: extracted['cat_no'] = parts[-1].replace('.html', '').replace('.htm', '')

            # Price Data
            price_data = await page.evaluate("""
                () => {
                    const rows = Array.from(document.querySelectorAll('table tr'));
                    const prices = [];
                    rows.forEach(row => {
                        const cols = row.querySelectorAll('td');
                        if (cols.length >= 3) {
                            const size = cols[0].innerText.trim();
                            const price = cols[1].innerText.trim();
                            if ((price.includes('$') || price.includes('€')) && !price.includes('Login')) {
                                prices.push({ size, price, availability: cols[2].innerText.trim() });
                            }
                        }
                    });
                    return prices;
                }
            """)

            return {
                "ambeed_cat_no": extracted['cat_no'],
                "cas_number": extracted['cas_no'],
                "product_name": title.split("|")[0].strip(),
                "product_url": url,
                "category": category,
                "smiles_code": extracted['smiles'],
                "properties": {
                    "formula": extracted['formula'],
                    "mw": extracted['mw'],
                    "purity": extracted['purity']
                },
                "price_data": price_data,
                "body_text": body_text,
                "source_name": "Ambeed",
                "crawled_at": datetime.utcnow().isoformat()
            }
        except Exception as e:
            logger.error(f"Extraction logic failed: {e}")
            return None

    async def _fetch_pubchem_smiles(self, cas_number: str) -> Optional[str]:
        """Fetch SMILES from PubChem API using CAS number"""
        if not cas_number: return None
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/property/IsomericSMILES/JSON"
            async with aiohttp.ClientSession() as session:
                async with session.get(url, timeout=10) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data['PropertyTable']['Properties'][0]['IsomericSMILES']
        except:
            return None

    async def _enrich_with_gemini(self, description: str, current_data: Dict) -> Dict:
        """Extract Target and Attributes using Gemini"""
        # Quick check to skip expensive calls
        if all([current_data.get('cas_number'), current_data.get('smiles_code'), current_data.get('properties', {}).get('mw')]):
             if len(description) < 200:
                 return {"target": None, "properties": {}, "summary": "Structured data extracted."}

        prompt = f"""
        Extract ADC-relevant data from this chemical reagent text.
        Text: {description[:2500]}
        
        Output JSON:
        {{
            "target": "Target Antigen (e.g. HER2) or null",
            "properties": {{ "purity": "...", "mw": "..." }},
            "summary": "1-sentence summary"
        }}
        """
        try:
            response = await self.model.generate_content_async(
                prompt, 
                generation_config=genai.GenerationConfig(response_mime_type="application/json")
            )
            return json.loads(response.text)
        except:
            return {}

    async def _get_embedding(self, text: str) -> Optional[List[float]]:
        try:
            return await rag_service.generate_embedding(text)
        except:
            return None

    async def _process_batch_enrichment(self, batch_data: List[Dict]):
        """Process a batch of items with AI in parallel"""
        if not batch_data: return
        tasks = [self._enrich_and_save_single(item) for item in batch_data]
        await asyncio.gather(*tasks)

    async def _enrich_and_save_single(self, raw_data: Dict):
        """Enrich a single item and save to DB"""
        try:
            body_text = raw_data.pop("body_text", "")
            
            # 1. PubChem Fallback
            if not raw_data.get('smiles_code') and raw_data.get('cas_number'):
                raw_data['smiles_code'] = await self._fetch_pubchem_smiles(raw_data['cas_number'])

            # 2. AI Enrichment
            ai_data = await self._enrich_with_gemini(body_text, raw_data)
            
            raw_data["target"] = ai_data.get("target")
            raw_data["summary"] = ai_data.get("summary")
            
            # Merge props
            if ai_data.get("properties"):
                raw_data["properties"].update(ai_data["properties"])

            # 3. Embedding
            embed_text = f"{raw_data['product_name']} {raw_data.get('target', '')} {raw_data.get('smiles_code', '')}"
            raw_data["embedding"] = await self._get_embedding(embed_text)

            # 4. Save
            res = supabase.table("commercial_reagents").upsert(raw_data, on_conflict="ambeed_cat_no").execute()
            
            if res.data:
                # Fire & Forget Refinement
                asyncio.create_task(self._trigger_refinement(res.data[0]))
                
        except Exception as e:
            logger.error(f"Save failed for {raw_data.get('product_name')}: {e}")

    async def _trigger_refinement(self, record: Dict):
        try:
            analysis = await ai_refiner.refine_single_record(record)
            if analysis and "error" not in analysis:
                update_data = {
                    "target": analysis.get("target"),
                    "ai_refined": True,
                    "properties": {**record.get("properties", {}), "ai_analysis": analysis}
                }
                supabase.table("commercial_reagents").update(update_data).eq("id", record["id"]).execute()
        except:
            pass

    async def run(self, search_term: str, limit: int, job_id: str):
        """Orchestrate the crawling process with Category Concurrency"""
        from app.api.scheduler import update_job_status
        await update_job_status(job_id, status="running")
        
        try:
            # Determine targets
            targets = {}
            if search_term and search_term != 'all':
                for cat, url in self.CATEGORIES.items():
                    if search_term.lower() in cat.lower():
                        targets[cat] = url
            else:
                targets = self.CATEGORIES

            if not targets:
                await update_job_status(job_id, status="completed", message="No categories matched.")
                return

            logger.info(f"🏁 Starting Concurrent Crawl for {len(targets)} categories...")
            
            # Create tasks for each category
            # Use semaphore to limit how many categories run strictly in parallel if needed, 
            # though we handled tab concurrency with global_semaphore.
            tasks = [self.crawl_category(cat, url, limit) for cat, url in targets.items()]
            
            results = await asyncio.gather(*tasks)
            total_processed = sum(results)
            
            await update_job_status(job_id, status="completed", records_drafted=total_processed, completed_at=datetime.utcnow().isoformat())

        except Exception as e:
            logger.error(f"Crawler run failed: {e}")
            await update_job_status(job_id, status="failed", errors=[str(e)])

# Singleton
ambeed_crawler = AmbeedCrawler()
