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
logger = logging.getLogger("CB_Antibody_Crawler")
logger.setLevel(logging.INFO)
if logger.hasHandlers(): logger.handlers.clear()

class FlushHandler(logging.FileHandler):
    def emit(self, record):
        super().emit(record)
        self.flush()

fh = FlushHandler("cb_antibody_crawl.log", encoding='utf-8')
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

class CB_AntibodyCrawler:
    BASE_URL = "https://www.creative-biolabs.com/bsab/"
    TARGETS_ROOT = "https://www.creative-biolabs.com/bsab/category/by-targets-1377.htm"

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

    def _rotate_proxy_index(self):
        self.current_proxy_index = (self.current_proxy_index + 1) % len(self.proxies)
        logger.info(f"🔄 Rotating Proxy to Index {self.current_proxy_index} (Port: {self.proxies[self.current_proxy_index]['server'].split(':')[-1]})")

    async def _init_browser(self, p) -> (Browser, BrowserContext):
        logger.info(f"🌐 Launching Browser with Local IP (Proxy Disabled)")
        
        browser = await p.chromium.launch(
            headless=True, # Headless for stability in long runs
            args=["--disable-blink-features=AutomationControlled"]
        )

        context = await browser.new_context(
            user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            viewport={"width": 1280, "height": 720}
        )
        return browser, context

    async def run(self, limit: int = 10000):
        logger.info(f"🚀 [CB ANTIBODY] Starting FULL Harvest (A-Z) via Local IP. Limit: {limit}")
        
        async with async_playwright() as p:
            browser, context = await self._init_browser(p)
            
            # 1. Fetch ALL Target Links from Root Page
            logger.info(f"🔗 Fetching Target List from {self.TARGETS_ROOT}")
            all_target_links = []
            
            try:
                page = await context.new_page()
                await page.goto(self.TARGETS_ROOT, timeout=60000)
                await asyncio.sleep(5)
                
                # Parse all links under 'hot-target' lists (A, B, C...)
                links = await page.evaluate("""
                    () => {
                        const anchors = Array.from(document.querySelectorAll('.hot-target li a'));
                        return anchors.map(a => a.href).filter(h => h.includes('symbolsearch-'));
                    }
                """)
                all_target_links = sorted(list(set(links)))
                logger.info(f"✅ Found {len(all_target_links)} unique target categories.")
                await page.close()
            except Exception as e:
                logger.error(f"❌ Failed to fetch target list: {e}")
                return

            total_processed = 0
            
            # 2. Iterate All Targets
            for t_url in all_target_links:
                if total_processed >= limit: break
                
                target_name = t_url.split('symbolsearch-')[-1].replace('.htm', '')
                logger.info(f"💓 [HEARTBEAT] Scanning Target: {target_name}")
                
                try:
                    page = await context.new_page()
                    await page.goto(t_url, timeout=60000)
                    await asyncio.sleep(random.uniform(5, 10)) # Increased delay
                    
                    product_links = await page.evaluate("""
                        () => {
                            const links = Array.from(document.querySelectorAll('a'));
                            return links
                                .map(a => a.href)
                                .filter(href => href.includes('bispecific-') && href.endsWith('.htm') && !href.includes('category/') && /\-\d+\.htm$/.test(href));
                        }
                    """)
                    product_links = list(set(product_links))
                    logger.info(f"🔎 Found {len(product_links)} products for {target_name}")
                    await page.close()
                    
                    # 3. Process Products
                    for p_url in product_links:
                        if total_processed >= limit: break
                        
                        success = await self.process_product(context, p_url)
                        
                        if success:
                            total_processed += 1
                            logger.info(f"📦 Total Saved: {total_processed}")
                        
                        # High Safety Delay for Local IP
                        delay = random.uniform(10, 20)
                        logger.info(f"💤 Sleeping {delay:.1f}s (Local IP Safety)...")
                        await asyncio.sleep(delay)
                        
                except Exception as e:
                    logger.error(f"Error scanning target {t_url}: {e}")
                    await asyncio.sleep(30) # Error penalty delay

            if browser: await browser.close()

    async def process_product(self, context, url):
        # 중복 체크
        cat_no_raw = url.split('-')[-1].replace('.htm', '')
        
        # Validate that cat_no_raw is a number (product ID)
        if not cat_no_raw.isdigit():
            logger.warning(f"⚠️ Invalid Product ID in URL: {url} (ID: {cat_no_raw}). Skipping.")
            return False

        # Handle cases where cat_no might be extracted differently
        # Usually CB-xxxx. Let's just use the url end part as key.
        cat_no = f"CB-{cat_no_raw}"
        
        existing = supabase.table("antibody_library").select("id").eq("cat_no", cat_no).execute()
        if existing.data: 
            logger.info(f"⏭️ Skipping existing: {cat_no}")
            return True

        page = await context.new_page()
        try:
            await page.goto(url, timeout=60000)
            # await asyncio.sleep(random.uniform(1, 3)) # Short wait for DOM
            
            title = await page.title()
            product_name = title.split('|')[0].strip()
            
            # --- NEW PARSING LOGIC (DIV BASED) ---
            raw_data = await page.evaluate("""
                () => {
                    const specs = {};
                    const targets_data = {
                        gene_ids: [],
                        uniprot_ids: []
                    };
                    
                    // 1. Summary Extraction
                    const summaryEl = document.querySelector('.showlistline .col-md-10');
                    const summary = summaryEl ? summaryEl.innerText.trim() : "";

                    // 2. Specs & Targets Extraction (Iterate over div-based rows)
                    const rows = document.querySelectorAll('.table-row-group');
                    
                    rows.forEach(row => {
                        const keyEl = row.querySelector('.table-cell-left');
                        const valEl = row.querySelector('.table-cell');
                        
                        if (keyEl && valEl) {
                            const key = keyEl.innerText.trim().replace(/:$/, '');
                            const val = valEl.innerText.trim();
                            
                            specs[key] = val;

                            // Extract Gene IDs
                            if (key === 'Gene ID') {
                                const links = valEl.querySelectorAll('a');
                                links.forEach(a => {
                                    const id = a.innerText.trim();
                                    if(id) targets_data.gene_ids.push(id);
                                });
                                if (links.length === 0 && val) {
                                    val.split(/[,;]/).forEach(t => targets_data.gene_ids.push(t.trim()));
                                }
                            }

                            // Extract UniProt IDs
                            if (key === 'UniProt ID') {
                                const links = valEl.querySelectorAll('a');
                                links.forEach(a => {
                                    const id = a.innerText.trim();
                                    if(id) targets_data.uniprot_ids.push(id);
                                });
                                if (links.length === 0 && val) {
                                    val.split(/[,;]/).forEach(t => targets_data.uniprot_ids.push(t.trim()));
                                }
                            }
                        }
                    });
                    
                    // 3. Construct Targets List
                    const targets = [];
                    const targetStr = specs['Targets'] || "";
                    const symbols = targetStr.split('&').map(s => s.trim()).filter(s => s);
                    
                    symbols.forEach((sym, index) => {
                        targets.push({
                            symbol: sym,
                            gene_id: targets_data.gene_ids[index] || null, 
                            uniprot_id: targets_data.uniprot_ids[index] || null,
                            role: `Target ${index + 1}`
                        });
                    });

                    return { specs, targets, summary };
                }
            """)
            
            # DB Logic
            ab_data = {
                "product_name": product_name,
                "cat_no": cat_no,
                "antibody_format": raw_data["specs"].get("Type"),
                "host_species": raw_data["specs"].get("Host Animal 1") or raw_data["specs"].get("Host"),
                "isotype": raw_data["specs"].get("Isotype"),
                "related_disease": raw_data["specs"].get("Related Disease"),
                "full_spec": raw_data["specs"],
                "source_url": url,
                "summary": raw_data["summary"] or f"Antibody targeting {raw_data['specs'].get('Targets')}"
            }
            
            # AI Embedding (Optional)
            if self.model and ab_data["summary"]:
                try:
                    embed_res = genai.embed_content(model="models/embedding-001", content=f"{product_name} {ab_data['summary']}")
                    ab_data["embedding"] = embed_res['embedding']
                except: pass

            # --- UPSERT WITH FALLBACK ---
            ab_res = supabase.table("antibody_library").upsert(ab_data, on_conflict="cat_no").execute()
            if not ab_res.data:
                # Fallback fetch
                ab_res = supabase.table("antibody_library").select("id").eq("cat_no", cat_no).execute()
            
            if not ab_res.data:
                logger.error(f"❌ DB Save Failed for {cat_no}")
                return False
            ab_id = ab_res.data[0]['id']
            
            # Targets
            for i, t_info in enumerate(raw_data["targets"]):
                symbol = t_info.get("symbol")
                if not symbol: continue
                
                t_master_data = {
                    "target_symbol": symbol,
                    "gene_id": int(t_info.get("gene_id")) if t_info.get("gene_id") and t_info.get("gene_id").isdigit() else None,
                    "uniprot_id": t_info.get("uniprot_id"),
                    "alternative_symbols": t_info.get("alt_names")
                }
                
                t_res = supabase.table("target_master").upsert(t_master_data, on_conflict="target_symbol").execute()
                if not t_res.data:
                    t_res = supabase.table("target_master").select("id").eq("target_symbol", symbol).execute()
                
                if t_res.data:
                    t_id = t_res.data[0]['id']
                    supabase.table("application_map").insert({
                        "antibody_id": ab_id,
                        "target_id": t_id,
                        "target_role": t_info.get("role"),
                        "interaction_type": "Primary Binder"
                    }).execute()

            logger.info(f"✅ Saved: {cat_no}")
            return True

        except Exception as e:
            logger.error(f"Error processing {url}: {e}")
            return False
        finally:
            await page.close()

if __name__ == "__main__":
    import sys
    # Default high limit for full crawl
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 10000 
    crawler = CB_AntibodyCrawler()
    asyncio.run(crawler.run(limit))
