import asyncio
import random
import logging
import json
import time
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

# --- Crawler Implementation ---

class AmbeedCrawlerLite:
    BASE_URL = "https://www.ambeed.com"
    
    # Categories with Max Pages (Approximation or safe high number)
    CATEGORIES = {
        "ADC Linker": {"url": "https://www.ambeed.com/adc-linkers.html", "max_page": 25},
        "ADC Toxins": {"url": "https://www.ambeed.com/adc-toxins.html", "max_page": 15},
        "ADC Payload-Linker": {"url": "https://www.ambeed.com/payload-linker-conjugate.html", "max_page": 20}
    }
    
    def __init__(self):
        self.ua = ua
        self.global_semaphore = asyncio.Semaphore(5) # Conservative concurrency

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
        return len(smiles) > 5

    async def _exists_in_db(self, cat_no: str) -> bool:
        """Check if catalog number already exists in DB"""
        try:
            res = supabase.table("commercial_reagents").select("id").eq("ambeed_cat_no", cat_no).limit(1).execute()
            return bool(res.data)
        except:
            return False

    async def process_page_items(self, context, page_url, category_name):
        """Process a single list page and return items"""
        page = await context.new_page()
        items = []
        try:
            # Random delay before page load
            await asyncio.sleep(random.uniform(2, 4))
            
            logger.info(f"📄 Visiting List: {page_url}")
            await page.goto(page_url, wait_until="domcontentloaded", timeout=60000)
            await asyncio.sleep(1) # Render wait

            # Extract Product Links
            products = await page.evaluate("""
                () => {
                    const links = Array.from(document.querySelectorAll('a[href*="/products/"]'));
                    return links
                        .map(a => ({ href: a.href, text: a.innerText }))
                        .filter(x => x.href.endsWith('.html') && !x.href.includes('category'));
                }
            """)
            
            # Deduplicate by URL
            unique_products = {p['href']: p for p in products}.values()
            
            for prod in unique_products:
                item_url = prod['href']
                cat_no_match = re.search(r'/products/([a-zA-Z0-9-]+)\.html', item_url)
                cat_no_candidate = cat_no_match.group(1) if cat_no_match else item_url.split('/')[-1].replace('.html', '')
                
                # DB Check
                if await self._exists_in_db(cat_no_candidate):
                    logger.info(f"⏩ [SKIP] {cat_no_candidate} exists.")
                    continue
                
                # Detail Page Processing
                details = await self._process_single_product(context, item_url, category_name)
                if details:
                    items.append(details)
                    # Short sleep between items to be polite
                    await asyncio.sleep(random.uniform(1, 2))
                    
        except Exception as e:
            logger.error(f"❌ Error processing list page {page_url}: {e}")
        finally:
            await page.close()
            
        return items, len(unique_products)

    async def _process_single_product(self, context, url, category):
        page = await context.new_page()
        try:
            await page.goto(url, wait_until="domcontentloaded", timeout=45000)
            
            # Parsing Logic based on User's Description
            # Expecting tables with CAS No., Formula, SMILES Code, MDL No.
            extracted = await page.evaluate("""
                () => {
                    const data = {
                        cas_no: null,
                        formula: null,
                        mw: null,
                        smiles: null,
                        mdl_no: null,
                        cat_no: null
                    };
                    
                    // Helper to clean text
                    const clean = (text) => text ? text.replace(/:/g, '').trim() : null;

                    // Strategy 1: Look for table rows (tr > td:key, td:value)
                    const rows = document.querySelectorAll('tr');
                    rows.forEach(row => {
                        const cells = row.querySelectorAll('td');
                        if (cells.length >= 2) {
                            const key = cells[0].innerText.trim();
                            const val = cells[1].innerText.trim();
                            
                            if (key.includes('CAS No')) data.cas_no = clean(val);
                            if (key.includes('Formula')) data.formula = clean(val);
                            if (key.includes('M.W')) data.mw = clean(val);
                            if (key.includes('SMILES Code')) data.smiles = clean(val);
                            if (key.includes('MDL No')) data.mdl_no = clean(val);
                        }
                    });

                    // Strategy 2: Fallback to scanning all elements if table fails (rare)
                    if (!data.cas_no) {
                        const allText = document.body.innerText;
                        const casMatch = allText.match(/CAS No\.?\s*:?\s*([0-9-]+)/i);
                        if (casMatch) data.cas_no = casMatch[1];
                    }
                    
                    // Get Product Name
                    const title = document.title.split('|')[0].trim();
                    const nameMatch = title.match(/Product Details of\s*\[?\s*([^\]]+)\s*\]?/i) || [null, title];
                    
                    return { ...data, product_name: nameMatch[1] || title };
                }
            """)
            
            # Cat No Extraction from URL if not found
            cat_no = url.split('/')[-1].replace('.html', '')
            
            final_data = {
                "ambeed_cat_no": cat_no,
                "product_name": extracted['product_name'],
                "product_url": url,
                "category": category,
                "cas_number": extracted['cas_no'],
                "formula": extracted['formula'],
                "molecular_weight": extracted['mw'],
                "smiles_code": extracted['smiles'],
                # "mdl_number": extracted['mdl_no'], # Removing top-level to avoid Schema Error
                "source_name": "Ambeed",
                "crawled_at": datetime.utcnow().isoformat(),
                "ai_refined": False,
                "properties": {
                    "mdl_number": extracted['mdl_no']
                }
            }
            
            logger.info(f"✅ Extracted: {cat_no} | {final_data['product_name']}")
            return final_data

        except Exception as e:
            logger.error(f"❌ Failed to process details {url}: {e}")
            return None
        finally:
            await page.close()

    async def _save_batch(self, items: List[Dict]):
        if not items: return
        try:
            logger.info(f"📤 Saving {len(items)} items to DB...")
            res = supabase.table("commercial_reagents").upsert(items, on_conflict="ambeed_cat_no").execute()
            if not res.data:
                 logger.warning("⚠️ Saved but no data returned (RLS?)")
        except Exception as e:
            logger.error(f"❌ DB Save Failed: {e}")

    async def run_round_robin(self):
        """
        Round Robin Strategy:
        Iterate through categories, fetching one page at a time.
        Linker Pg1 -> Toxin Pg1 -> Conjugate Pg1 -> Linker Pg2 ...
        """
        
        # State tracking for each category
        category_states = {
            name: {"page": 1, "done": False, "url": meta["url"], "max": meta["max_page"], "retry_count": 0}
            for name, meta in self.CATEGORIES.items()
        }
        
        async with async_playwright() as p:
            context = await self._init_browser(p)
            
            round_num = 1
            while True:
                active_cats = [n for n, s in category_states.items() if not s["done"]]
                if not active_cats:
                    logger.info("🎉 All categories completed!")
                    break
                
                logger.info(f"=== 🔄 Starting Round {round_num} (Active: {active_cats}) ===")
                
                for cat_name in active_cats:
                    state = category_states[cat_name]
                    current_page = state["page"]
                    
                    if current_page > state["max"]:
                        logger.info(f"🛑 {cat_name} reached max page limit ({state['max']}). Marking done.")
                        state["done"] = True
                        continue
                        
                    # Construct URL
                    base_url = state["url"]
                    separator = "&" if "?" in base_url else "?"
                    page_url = f"{base_url}{separator}pagesize=20&pageindex={current_page}"
                    
                    logger.info(f"🎯 [{cat_name}] Fetching Page {current_page} (Attempt {state['retry_count'] + 1}/3)...")
                    
                    # Fetch Items (Modified to unpack tuple)
                    items, raw_count = await self.process_page_items(context, page_url, cat_name)
                    
                    if raw_count == 0:
                        # Possible block or temporary empty page
                        state["retry_count"] += 1
                        logger.warning(f"⚠️ {cat_name} Page {current_page} returned 0 items. Retry count: {state['retry_count']}/3")
                        
                        if state["retry_count"] >= 3:
                            logger.warning(f"🛑 {cat_name} failed 3 times consecutively. Skipping Page {current_page} and moving to next.")
                            state["page"] += 1
                            state["retry_count"] = 0
                        else:
                            logger.info(f"⏳ Keeping {cat_name} active to retry Page {current_page} in next round.")
                            # Do NOT increment page, so we retry same page next time
                    
                    elif not items and raw_count > 0:
                        # Found items but all were duplicates
                        logger.info(f"⏭️ {cat_name} Page {current_page}: All {raw_count} items exist in DB. Moving to next page.")
                        state["page"] += 1
                        state["retry_count"] = 0 # Reset retry on valid page load
                    else:
                        # Found new items
                        await self._save_batch(items)
                        state["page"] += 1
                        state["retry_count"] = 0 # Reset retry on success
                    
                    # Inter-category delay to avoid blocking
                    await asyncio.sleep(random.uniform(5, 10))
                
                round_num += 1
                # Inter-round delay
                logger.info("💤 Resting between rounds...")
                await asyncio.sleep(10)

if __name__ == "__main__":
    crawler = AmbeedCrawlerLite()
    asyncio.run(crawler.run_round_robin())
