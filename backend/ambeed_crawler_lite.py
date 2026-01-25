import asyncio
import random
import logging
import json
import time
import subprocess
import os
import re
from datetime import datetime
from typing import Dict, List, Optional, Any
from dotenv import load_dotenv

# Load env immediately
load_dotenv()

# Setup Logging
logger = logging.getLogger("AmbeedCrawlerLite")
logger.setLevel(logging.INFO)

# Clear existing handlers
if logger.hasHandlers():
    logger.handlers.clear()

class FlushHandler(logging.FileHandler):
    def emit(self, record):
        super().emit(record)
        self.flush()

# File Handler with explicit flush
fh = FlushHandler("ambeed_crawl.log", encoding='utf-8')
fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(fh)

# Console Handler
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(ch)

# Debug Env
url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_KEY")
if not url or not key:
    logger.error("❌ SUPABASE_URL or SUPABASE_SERVICE_KEY missing in environment variables!")
else:
    logger.info(f"✅ Supabase Env Found: URL={url}, KEY={'*' * 10}")

# --- Imports with Fallbacks ---

try:
    from playwright.async_api import async_playwright, Page, BrowserContext, TimeoutError as PlaywrightTimeoutError
except ImportError:
    logger.error("Playwright is missing. Please install it.")
    raise

try:
    from fake_useragent import UserAgent
    ua = UserAgent()
except ImportError:
    class MockUA:
        random = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    ua = MockUA()

try:
    import google.generativeai as genai
    HAS_GENAI = True
except ImportError:
    logger.warning("google.generativeai missing.")
    HAS_GENAI = False

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    logger.warning("rdkit missing. SMILES validation disabled.")
    HAS_RDKIT = False

try:
    import aiohttp
except ImportError:
    pass

# Mock Heavy Dependencies
try:
    from app.core.config import settings
    from app.core.supabase import supabase
    from supabase import create_client, Client # Import explicitly for re-init
    
    # Check if we got a mock client and try to fix it
    if hasattr(supabase, "is_mock") or (hasattr(supabase, "execute") and "MockSupabaseClient" in str(supabase.execute)):
        logger.warning("⚠️ Imported Supabase client is Mock. Attempting to re-initialize with env vars...")
        _url = os.getenv("SUPABASE_URL")
        _key = os.getenv("SUPABASE_SERVICE_KEY")
        if _url and _key:
            try:
                supabase = create_client(_url, _key)
                logger.info("✅ Supabase client re-initialized successfully.")
            except Exception as e:
                logger.error(f"❌ Failed to re-initialize Supabase: {e}")

except ImportError:
    class MockSettings:
        GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY", "")
    settings = MockSettings()
    
    class MockSupabase:
        def table(self, name): return self
        def upsert(self, *args, **kwargs): return self
        def execute(self): return type('obj', (object,), {'data': []})
    
    # Try to init real client if imports allow
    try:
        from supabase import create_client
        _url = os.getenv("SUPABASE_URL")
        _key = os.getenv("SUPABASE_SERVICE_KEY")
        if _url and _key:
             supabase = create_client(_url, _key)
        else:
             supabase = MockSupabase()
    except:
        supabase = MockSupabase()

# Mock RAG Service
class MockRagService:
    async def generate_embedding(self, text):
        return [0.1] * 768 # Dummy embedding
rag_service = MockRagService()

# Mock Scheduler
async def update_job_status(*args, **kwargs): pass
async def is_cancelled(*args, **kwargs): return False
async def get_job_from_db(*args, **kwargs): return {}


# --- Crawler Implementation ---

class AmbeedCrawlerLite:
    BASE_URL = "https://www.ambeed.com"
    CATEGORIES = {
        "ADC Toxins": "https://www.ambeed.com/adc-toxins.html",
        "ADC Linker": "https://www.ambeed.com/adc-linkers.html",
        "ADC Payload-Linker": "https://www.ambeed.com/payload-linker-conjugate.html"
    }
    
    def __init__(self):
        self.ua = ua
        self.global_semaphore = asyncio.Semaphore(1)
        self.model = None

    def _get_model(self):
        if not HAS_GENAI: return None
        if not self.model:
            if not settings.GOOGLE_API_KEY:
                logger.warning("GOOGLE_API_KEY not found in settings.")
                return None
            genai.configure(api_key=settings.GOOGLE_API_KEY)
            self.model = genai.GenerativeModel('gemini-2.0-flash')
        return self.model

    async def _init_browser(self, p) -> BrowserContext:
        try:
            browser = await p.chromium.launch(
                headless=True,
                args=['--no-sandbox', '--disable-setuid-sandbox', '--disable-dev-shm-usage', '--disable-gpu']
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
            logger.error(f"Browser Launch Failed: {e}")
            raise e

    def validate_smiles(self, smiles: str) -> bool:
        if not smiles: return False
        if HAS_RDKIT:
            try:
                mol = Chem.MolFromSmiles(smiles)
                return mol is not None
            except: 
                return False
        return len(smiles) > 5 # Simple length check if RDKit missing

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10, job_id: str = None, start_page: int = 1, batch_size: int = 2) -> int:
        logger.info(f"🚀 [AMBEED LITE] {category_name} (Start Page: {start_page}, Limit: {limit})")
        
        count = 0
        page_num = start_page
        batch_data = []
        seen_hrefs = set() # Track duplicates in this session
        
        async with async_playwright() as p:
            try:
                context = await self._init_browser(p)
                
                while count < limit:
                    # 1. Clear Cookies to appear as a new guest
                    await context.clear_cookies()
                    page = await context.new_page()
                    
                    # 2. Random Delay (1.5s ~ 3.5s)
                    delay = random.uniform(1.5, 3.5)
                    await asyncio.sleep(delay)
                    
                    separator = "&" if "?" in base_url else "?"
                    url = base_url if page_num == 1 else f"{base_url}{separator}pagesize=20&pageindex={page_num}"
                    logger.info(f"💓 [CRAWLER_HEARTBEAT] {datetime.now().isoformat()} - Category: {category_name}, Page: {page_num}")
                    logger.info(f"🌐 [PAGE {page_num}] Accessing: {url} (Sleep: {delay:.2f}s)")
                    
                    try:
                        await page.goto(url, wait_until="domcontentloaded", timeout=60000)
                        await asyncio.sleep(1) # Extra buffer for JS rendering
                        
                        products = await page.evaluate("""
                            () => {
                                let items = Array.from(document.querySelectorAll('.product-item, .item, .product-info-main, .product-item-info'));
                                if (items.length === 0) {
                                    return Array.from(document.querySelectorAll('a[href*="/products/"]'))
                                        .map(a => ({ href: a.href, cat_no: null }))
                                        .filter((v, i, a) => a.findIndex(t => t.href === v.href) === i);
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
                            links = await page.evaluate("""
                                () => Array.from(document.querySelectorAll('a[href*="/products/"]'))
                                    .map(a => a.href)
                                    .filter(href => href.includes('.html'))
                            """)
                            products = [{"href": link, "cat_no": None} for link in set(links)]

                        if not products:
                            logger.warning(f"⚠️ [PAGE {page_num}] No products found. Likely end of category.")
                            break
                        
                        # 3. Duplicate Detection and Logging
                        new_items = [p for p in products if p["href"] not in seen_hrefs]
                        dup_count = len(products) - len(new_items)
                        
                        if dup_count > 0:
                            logger.info(f"⏭️ [SKIP] {dup_count} items are duplicates. (Unique in this page: {len(new_items)})")
                        
                        if not new_items:
                            logger.warning(f"🚩 [STUCK] Page {page_num} returned ONLY duplicates. Forcing next page...")
                        else:
                            logger.info(f"📦 [PAGE {page_num}] Processing {len(new_items)} new products.")

                        for prod in new_items:
                            if count >= limit: break
                            
                            link = prod["href"]
                            seen_hrefs.add(link)
                            cat_no = prod["cat_no"] or link.split('/')[-1].replace('.html', '')
                            
                            logger.info(f"🔨 Refining record: [Cat No.{cat_no}] (AI Enriching...)")

                            res = await self._process_single_product(context, link, category_name)
                            if res:
                                final_item = await self._enrich_and_prepare_item(res)
                                if final_item:
                                    batch_data.append(final_item)
                                    count += 1
                                    
                                    if len(batch_data) >= batch_size:
                                        await self._save_batch(batch_data)
                                        batch_data = []

                        # Force Jump to Next Page
                        logger.info(f"✅ [PAGE {page_num}] Completed. [CRAWLER] Moving to Page {page_num + 1}...")
                        page_num += 1
                        await page.close() # Close page to save memory

                    except Exception as e:
                        logger.error(f"❌ Error on page {page_num}: {e}")
                        logger.info(f"⚠️ [Recovery] Ignoring error. [CRAWLER] Moving to Page {page_num + 1}...")
                        page_num += 1
                        await page.close()
                        await asyncio.sleep(2)
                
                if batch_data:
                    await self._save_batch(batch_data)
            
            finally: 
                await context.close()
        return count
                
                if batch_data:
                    await self._save_batch(batch_data)
            
            finally: 
                await context.close()
        return count

    async def _process_single_product(self, context, url, category):
        async with self.global_semaphore:
            page = await context.new_page()
            try:
                await page.goto(url, wait_until="domcontentloaded", timeout=30000)
                title = await page.title()
                
                extracted = await page.evaluate("""
                    () => {
                        const results = { cat_no: '', cas_no: '', smiles: '' };
                        const tds = Array.from(document.querySelectorAll('td, th, span, div'));
                        for (let i = 0; i < tds.length; i++) {
                            const text = tds[i].innerText;
                            if (text.includes('SMILES')) {
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

    async def _enrich_and_prepare_item(self, raw_data):
        smiles = raw_data.get("smiles_code")
        is_valid = self.validate_smiles(smiles)
        
        # Simple AI enrichment if available
        ai_data = await self._enrich_with_gemini(raw_data, smiles)
        
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
            "payload_smiles": ai_data.get("payload_smiles"),
            "linker_smiles": ai_data.get("linker_smiles"),
            "full_smiles": ai_data.get("full_smiles") or smiles,
            "cas_number": raw_data.get("cas_number"),
            "target": ai_data.get("target"),
            "summary": ai_data.get("summary"),
            "properties": extended_properties,
            "crawled_at": raw_data["crawled_at"],
            "ai_refined": True
        }
        
        try:
             # Embedding
            embed_text = f"{final_data['product_name']} {final_data.get('smiles_code') or ''}"
            embedding = await rag_service.generate_embedding(embed_text)
            final_data["embedding"] = embedding
        except: pass
        
        return final_data

    async def _enrich_with_gemini(self, raw_data: Dict, smiles: Optional[str] = None) -> Dict:
        if not HAS_GENAI: return {}
        
        prompt = f"""
        Extract detailed ADC reagent information:
        Title: {raw_data['product_name']}
        Input SMILES: {smiles or 'None'}

        Output JSON:
        {{
            "payload_smiles": null,
            "linker_smiles": null,
            "full_smiles": null,
            "target": null,
            "summary": "Extracted from Ambeed"
        }}
        """
        try:
            model = self._get_model()
            if not model: return {}
            response = await model.generate_content_async(prompt, generation_config=genai.GenerationConfig(response_mime_type="application/json"))
            data = json.loads(response.text)
            if isinstance(data, list):
                return data[0] if data else {}
            return data
        except Exception as e:
            logger.error(f"Gemini API Error: {e}")
            return {}

    async def _save_batch(self, items: List[Dict]):
        try:
            if not items: return
            logger.info(f"📤 Saving {len(items)} items to Supabase...")
            res = supabase.table("commercial_reagents").upsert(items, on_conflict="ambeed_cat_no").execute()
            
            # Check response data
            if res.data:
                logger.info(f"✅ Saved Successfully! Inserted/Updated {len(res.data)} rows.")
                # Print first item ID/CatNo to confirm
                logger.info(f"   Sample: {res.data[0].get('ambeed_cat_no')} - {res.data[0].get('product_name')}")
            else:
                logger.warning("⚠️ HTTP 200 but NO DATA returned. Likely RLS policy or permission issue.")
                
        except Exception as e:
            logger.error(f"❌ Save failed: {e}")
            if hasattr(e, 'code'): logger.error(f"   Error Code: {e.code}")
            if hasattr(e, 'details'): logger.error(f"   Error Details: {e.details}")

    async def run(self, search_term: str, limit: int, job_id: str, start_page: int = 1, batch_size: int = 2):
        targets = {cat: url for cat, url in self.CATEGORIES.items() if not search_term or search_term == 'all' or search_term.lower() in cat.lower()}
        for cat, url in targets.items():
            await self.crawl_category(cat, url, limit, job_id, start_page, batch_size)

if __name__ == "__main__":
    import sys
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    category = sys.argv[2] if len(sys.argv) > 2 else "ADC Toxins"
    start_page = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    
    print(f"🚀 Starting Bulk Crawl: {category}, Limit: {limit}, Start Page: {start_page}")
    crawler = AmbeedCrawlerLite()
    asyncio.run(crawler.run(category, limit, "bulk_job_lite", start_page=start_page, batch_size=10))
