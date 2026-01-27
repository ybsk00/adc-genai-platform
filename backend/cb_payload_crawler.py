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
logger = logging.getLogger("CB_Payload_Crawler")
logger.setLevel(logging.INFO)
if logger.hasHandlers(): logger.handlers.clear()

class FlushHandler(logging.FileHandler):
    def emit(self, record):
        super().emit(record)
        self.flush()

fh = FlushHandler("cb_payload_crawl.log", encoding='utf-8')
fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(ch)

# --- Imports ---
try:
    from playwright.async_api import async_playwright, Browser, BrowserContext
except ImportError:
    logger.error("Playwright is missing.")
    raise

try:
    import google.generativeai as genai
    HAS_GENAI = True
except ImportError:
    HAS_GENAI = False

try:
    from app.core.supabase import supabase
    from supabase import create_client
    # Re-init if mock
    if hasattr(supabase, "is_mock") or "Mock" in str(type(supabase)):
        _url, _key = os.getenv("SUPABASE_URL"), os.getenv("SUPABASE_SERVICE_KEY")
        if _url and _key: supabase = create_client(_url, _key)
except ImportError:
    logger.error("Supabase client not found.")
    raise

class CB_PayloadCrawler:
    BASE_URL = "https://www.creative-biolabs.com/adc/classify-adc-toxins-6.htm"
    TOTAL_PAGES = 8

    # Proxy Configuration
    PROXY_HOST = "gate.decodo.com"
    PROXY_USER = "spa2akr23a"
    PROXY_PASS = "zgS51qwDcLAyc88a+c"
    PROXY_PORTS = range(10001, 10011) # 10001 to 10010

    def __init__(self):
        self.semaphore = asyncio.Semaphore(1)
        self.model = None
        if HAS_GENAI and os.getenv("GOOGLE_API_KEY"):
            genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))
            self.model = genai.GenerativeModel('gemini-1.5-flash')
        
        # Prepare Proxy List
        self.proxies = [
            {
                "server": f"http://{self.PROXY_HOST}:{port}",
                "username": self.PROXY_USER,
                "password": self.PROXY_PASS
            }
            for port in self.PROXY_PORTS
        ]
        self.current_proxy_index = 0

    def _get_current_proxy(self):
        proxy = self.proxies[self.current_proxy_index]
        return proxy

    async def _init_browser(self, p) -> (Browser, BrowserContext):
        logger.info(f"🌐 Launching Browser (Local IP via VPN)")
        
        browser = await p.chromium.launch(
            headless=True,
            args=["--disable-blink-features=AutomationControlled"]
        )

        context = await browser.new_context(
            user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            viewport={"width": 1280, "height": 720}
        )
        return browser, context

    async def run(self, limit: int = 10000):
        logger.info(f"🚀 [CB PAYLOAD] Starting FULL Harvest (VPN Mode). Limit: {limit}")
        
        async with async_playwright() as p:
            browser, context = await self._init_browser(p)
            
            total_processed = 0
            
            for page_num in range(1, self.TOTAL_PAGES + 1):
                if total_processed >= limit: break
                
                list_url = f"{self.BASE_URL}?page={page_num}"
                logger.info(f"📄 Processing Page {page_num}: {list_url}")
                
                try:
                    page = await context.new_page()
                    
                    product_links = []
                    for attempt in range(3):
                        await page.goto(list_url, timeout=60000)
                        content = await page.content()
                        if "Server is busy" in content:
                            logger.warning(f"⚠️ Server busy on list page {page_num}. Waiting 120s...")
                            await asyncio.sleep(120)
                            continue
                        
                        product_links = await page.evaluate("""
                            () => {
                                const anchors = Array.from(document.querySelectorAll('a'));
                                return anchors
                                    .map(a => a.href)
                                    .filter(h => h.includes('/adc/') && /\-\d+\.htm/.test(h) && !h.includes('classify-') && !h.includes('category'));
                            }
                        """)
                        if product_links: break
                        await asyncio.sleep(10)

                    product_links = sorted(list(set(product_links)))
                    logger.info(f"🔎 Found {len(product_links)} products on page {page_num}")
                    await page.close()
                    
                    for p_url in product_links:
                        if total_processed >= limit: break
                        
                        success = await self.process_product(context, p_url)
                        if success:
                            total_processed += 1
                        
                        # High Safety Delay (10-20s)
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
                await page.goto(url, timeout=60000)
                content = await page.content()
                if "Server is busy" in content:
                    logger.warning(f"⚠️ Server busy on product. Waiting 60s...")
                    await asyncio.sleep(60)
                    continue
                break

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
                    if(!data.title) data.title = document.title; // Fallback

                    // Strategy 1: .proshow-unit ul (Standard)
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

                    return data;
                }
            """)
            
            specs = raw_data['specs']
            raw_title = raw_data['title']
            
            # Clean name
            if 'CAT#:' in raw_title:
                product_name = raw_title.split('CAT#:')[0]
            elif '(' in raw_title:
                 product_name = raw_title.split('(')[0]
            else:
                product_name = raw_title
                
            product_name = product_name.replace(" - Creative Biolabs", "").strip()
            
            # 1. CAT NO Extraction Strategy
            cat_no = specs.get('Product No') or specs.get('Cat. No') or specs.get('Catalog Number')
            
            if not cat_no:
                # Regex from title
                match = re.search(r'CAT#:\s*([A-Za-z0-9-]+)', raw_title, re.IGNORECASE)
                if match:
                    cat_no = match.group(1).strip()
            
            if not cat_no:
                match = re.search(r'-(\d+)\.htm$', url)
                if match:
                    cat_no = f"CBL-ADC-{match.group(1)}"
                    logger.warning(f"⚠️ Cat No not found. Generated: {cat_no}")
                else:
                    logger.error(f"❌ Could not determine Cat No for {url}")
                    return False
            
            # 2. CAS NO Extraction
            cas_number = None
            for k, v in specs.items():
                if k.lower() in ['cas no', 'cas number', 'cas #', 'cas']:
                    cas_number = v
                    break
            
            # If still missing CAS, try to find in specs keys loosely
            if not cas_number:
                # specific check for 'CAS NO' in keys
                pass

            classification = specs.get('Classification') or 'ADC Toxins'
            mw = specs.get('Molecular Weight') or specs.get('MW')
            
            properties = {}
            exclude_keys = ['product no', 'cat. no', 'catalog number', 'cas no', 'cas number', 'classification', 'molecular weight', 'mw']
            for k, v in specs.items():
                if k.lower() not in exclude_keys:
                    properties[k] = v
            
            record = {
                "ambeed_cat_no": cat_no,
                "product_name": product_name,
                "category": classification,
                "cas_number": cas_number,
                "molecular_weight": mw,
                "source_name": "Creative Biolabs",
                "product_url": url,
                "properties": properties,
                "crawled_at": datetime.now().isoformat(),
            }

            res = supabase.table("commercial_reagents").upsert(
                record, on_conflict="ambeed_cat_no"
            ).execute()
            
            if res.data or res:
                logger.info(f"✅ Saved: {cat_no} | CAS: {cas_number} | Name: {product_name}")
                return True

        except Exception as e:
            logger.error(f"Error processing {url}: {e}")
            return False
        finally:
            await page.close()

if __name__ == "__main__":
    crawler = CB_PayloadCrawler()
    asyncio.run(crawler.run())
