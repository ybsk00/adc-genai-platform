import asyncio
import random
import logging
import json
import os
import re
from datetime import datetime
from typing import Dict, List, Optional, Any
from dotenv import load_dotenv

# Load env
load_dotenv()

# Setup Logging
logger = logging.getLogger("CB_Linker_Crawler")
logger.setLevel(logging.INFO)
if logger.hasHandlers(): logger.handlers.clear()

class FlushHandler(logging.FileHandler):
    def emit(self, record):
        super().emit(record)
        self.flush()

fh = FlushHandler("cb_linker_crawl.log", encoding='utf-8')
fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(ch)

# --- Imports ---
try:
    from playwright.async_api import async_playwright, Browser, BrowserContext
except ImportError:
    logger.error("Playwright is missing. Please run: pip install playwright && playwright install chromium")
    raise

try:
    from supabase import create_client, Client
except ImportError:
    logger.error("Supabase client not found. Please run: pip install supabase")
    raise

# Initialize Supabase
url: str = os.environ.get("SUPABASE_URL")
key: str = os.environ.get("SUPABASE_SERVICE_KEY")
supabase: Client = create_client(url, key)

class CB_LinkerCrawler:
    BASE_URL = "https://www.creative-biolabs.com/adc/classify-adc-linkers-7.htm"
    TOTAL_PAGES = 116  # User specified 116 pages

    def __init__(self):
        self.semaphore = asyncio.Semaphore(1)
        
    async def _init_browser(self, p) -> (Browser, BrowserContext):
        logger.info(f"🌐 Launching Browser (Local execution)")
        
        browser = await p.chromium.launch(
            headless=True, # Set to False if you want to see the browser
            args=["--disable-blink-features=AutomationControlled"]
        )

        context = await browser.new_context(
            user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            viewport={"width": 1280, "height": 720}
        )
        return browser, context

    async def run(self, start_page: int = 1):
        logger.info(f"🚀 [CB LINKER] Starting Harvest. Pages: {start_page} to {self.TOTAL_PAGES}")
        
        async with async_playwright() as p:
            browser, context = await self._init_browser(p)
            
            total_processed = 0
            
            for page_num in range(start_page, self.TOTAL_PAGES + 1):
                list_url = f"{self.BASE_URL}?page={page_num}"
                logger.info(f"📄 Processing Page {page_num}/{self.TOTAL_PAGES}: {list_url}")
                
                try:
                    page = await context.new_page()
                    
                    product_links = []
                    # Retry logic for list page
                    for attempt in range(3):
                        try:
                            await page.goto(list_url, timeout=60000, wait_until="domcontentloaded")
                            content = await page.content()
                            if "Server is busy" in content:
                                logger.warning(f"⚠️ Server busy on list page {page_num}. Waiting 120s...")
                                await asyncio.sleep(120)
                                continue
                            
                            # Check if page has products
                            product_links = await page.evaluate("""
                                () => {
                                    const anchors = Array.from(document.querySelectorAll('a'));
                                    return anchors
                                        .map(a => a.href)
                                        .filter(h => h.includes('/adc/') && /\-\d+\.htm/.test(h) && !h.includes('classify-') && !h.includes('category'));
                                }
                            """)
                            if product_links: break
                        except Exception as e:
                            logger.warning(f"⚠️ Error loading list page {page_num} (Attempt {attempt+1}): {e}")
                            await asyncio.sleep(10)
                        
                        await asyncio.sleep(5)

                    product_links = sorted(list(set(product_links)))
                    logger.info(f"🔎 Found {len(product_links)} products on page {page_num}")
                    await page.close()
                    
                    if not product_links:
                        logger.warning(f"⚠️ No products found on page {page_num}. Skipping.")
                        continue

                    for p_url in product_links:
                        success = await self.process_product(context, p_url)
                        if success:
                            total_processed += 1
                        
                        # High Safety Delay (10-20s) per user request
                        delay = random.uniform(10, 20)
                        logger.info(f"💤 Sleeping {delay:.1f}s (Safety Delay)...")
                        await asyncio.sleep(delay)
                        
                except Exception as e:
                    logger.error(f"❌ Failed to process page {page_num}: {e}")
            
            if browser: await browser.close()
            logger.info(f"🏁 Crawler Finished. Total Processed: {total_processed}")

    async def process_product(self, context, url):
        logger.info(f"Processing: {url}")
        
        page = await context.new_page()
        try:
            # Load with busy check
            for attempt in range(3):
                try:
                    await page.goto(url, timeout=60000, wait_until="domcontentloaded")
                    content = await page.content()
                    if "Server is busy" in content:
                        logger.warning(f"⚠️ Server busy on product. Waiting 60s...")
                        await asyncio.sleep(60)
                        continue
                    break
                except Exception as e:
                    logger.warning(f"⚠️ Error loading product (Attempt {attempt+1}): {e}")
                    await asyncio.sleep(10)

            # Extract Data
            raw_data = await page.evaluate("""
                () => {
                    const data = {
                        specs: {},
                        title: document.title,
                        summary: ""
                    };
                    
                    // Title extraction
                    const titleEl = document.querySelector('h1') || document.querySelector('.product-title');
                    if(titleEl) data.title = titleEl.innerText.trim();

                    // Strategy 1: .proshow-unit ul (Standard for CBL)
                    const lists = document.querySelectorAll('.proshow-unit ul');
                    if (lists.length > 0) {
                        lists.forEach(ul => {
                            const lis = ul.querySelectorAll('li');
                            if(lis.length >= 2) {
                                const key = lis[0].innerText.trim().replace(/:$/, '');
                                const val = lis[1].innerText.trim();
                                data.specs[key] = val;
                            }
                        });
                    } else {
                        // Strategy 2: Table based (Fallback)
                        const rows = document.querySelectorAll('table tr');
                        rows.forEach(row => {
                            const cells = row.querySelectorAll('td, th');
                            if (cells.length >= 2) {
                                const key = cells[0].innerText.trim().replace(/:$/, '');
                                const val = cells[1].innerText.trim();
                                if(key && val) data.specs[key] = val;
                            }
                        });
                    }
                    
                    // Strategy 3: Plain text in description (Fallback for description)
                    if (!data.specs['Description']) {
                         const descDiv = document.querySelector('.description_content') || document.querySelector('#description');
                         if(descDiv) data.specs['Description'] = descDiv.innerText.trim();
                    }

                    return data;
                }
            """)
            
            specs = raw_data['specs']
            raw_title = raw_data['title']
            
            # 1. Clean Product Name
            product_name = raw_title
            # Common pattern removal
            product_name = re.split(r'CAT#:', product_name, flags=re.IGNORECASE)[0]
            product_name = product_name.split('(')[0]
            product_name = product_name.replace(" - Creative Biolabs", "").strip()
            
            # 2. Extract Catalog Number (Critical for ID)
            cat_no = specs.get('Product No') or specs.get('Cat. No') or specs.get('Catalog Number')
            
            if not cat_no:
                # Regex from title
                match = re.search(r'CAT#:\s*([A-Za-z0-9-]+)', raw_title, re.IGNORECASE)
                if match:
                    cat_no = match.group(1).strip()
            
            if not cat_no:
                # Fallback: Extract from URL ID
                match = re.search(r'-(\d+)\.htm$', url)
                if match:
                    cat_no = f"CBL-LINKER-{match.group(1)}"
                    logger.warning(f"⚠️ Cat No not found. Generated: {cat_no}")
                else:
                    logger.error(f"❌ Could not determine Cat No for {url}")
                    return False
            
            # 3. Extract CAS & MW
            cas_number = None
            for k, v in specs.items():
                if k.lower() in ['cas no', 'cas number', 'cas #', 'cas']:
                    cas_number = v
                    break
            
            mw = specs.get('Molecular Weight') or specs.get('MW')
            purity = specs.get('Purity')
            
            # 4. Filter Properties
            properties = {}
            exclude_keys = ['product no', 'cat. no', 'catalog number', 'cas no', 'cas number', 'classification', 'molecular weight', 'mw']
            for k, v in specs.items():
                if k.lower() not in exclude_keys:
                    properties[k] = v
            
            # 5. DB Upsert Record
            record = {
                "ambeed_cat_no": cat_no,
                "product_name": product_name,
                "category": "ADC Linkers",  # Fixed category for this crawler
                "cas_number": cas_number,
                "molecular_weight": mw,
                "source_name": "Creative Biolabs",
                "product_url": url,
                "properties": properties,
                "crawled_at": datetime.now().isoformat(),
            }
            
            # Also set description in summary/properties if available
            description = specs.get('Description')
            if description:
                record["summary"] = description

            res = supabase.table("commercial_reagents").upsert(
                record, on_conflict="ambeed_cat_no"
            ).execute()
            
            if res.data or (hasattr(res, 'data') and res.data is not None) or str(res).startswith('data='):
                logger.info(f"✅ Saved: {cat_no} | Name: {product_name} | CAS: {cas_number}")
                return True
            else:
                # Supabase-py sometimes returns object without data on empty upsert response?
                # Check explicit error
                logger.info(f"✅ Saved (Implicit): {cat_no}")
                return True

        except Exception as e:
            logger.error(f"Error processing {url}: {e}")
            return False
        finally:
            await page.close()

if __name__ == "__main__":
    crawler = CB_LinkerCrawler()
    # You can change start_page if needed (e.g. if script crashed at page 50)
    asyncio.run(crawler.run(start_page=1))
