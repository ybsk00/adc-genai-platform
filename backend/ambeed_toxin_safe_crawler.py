import asyncio
import random
import logging
import json
import os
import re
from datetime import datetime
from typing import Dict, List, Optional
from dotenv import load_dotenv

# Load env
load_dotenv()

# Setup Logging
logger = logging.getLogger("Ambeed_Toxin_Safe")
logger.setLevel(logging.INFO)
if logger.hasHandlers(): logger.handlers.clear()

fh = logging.FileHandler("ambeed_toxin_safe.log", encoding='utf-8')
fh.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(fh)

ch = logging.StreamHandler()
ch.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logger.addHandler(ch)

# Imports
try:
    from playwright.async_api import async_playwright, BrowserContext
except ImportError:
    logger.error("Playwright is missing.")
    raise

try:
    from supabase import create_client, Client
    url = os.getenv("SUPABASE_URL")
    key = os.getenv("SUPABASE_SERVICE_KEY")
    supabase: Client = create_client(url, key)
except ImportError:
    logger.error("Supabase client missing.")
    raise

class AmbeedToxinSafeCrawler:
    # URL 변경: ADC Toxins
    BASE_URL = "https://www.ambeed.com/adc-toxins.html"
    
    def __init__(self):
        self.page_start = 1
        self.target_page = 6 # 5-6페이지까지 시도
        self.batch_size = 2
        self.rest_duration = 120
        self.max_pages_limit = 20 # 넉넉하게 설정하되 막히면 중단
        
    async def _init_browser(self, p) -> (any, BrowserContext):
        logger.info(f"🌐 Launching Browser (Headed Mode) for Toxins")
        browser = await p.chromium.launch(
            headless=False, # 눈에 보이게 실행
            slow_mo=200,
            args=[
                '--no-sandbox',
                '--disable-blink-features=AutomationControlled',
                '--window-size=1920,1080',
            ]
        )
        context = await browser.new_context(
            user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/121.0.0.0 Safari/537.36",
            viewport={'width': 1920, 'height': 1080}
        )
        return browser, context

    async def wait_for_user_input(self, prompt="Press Enter to continue..."):
        print(f"\n👉 {prompt}")
        await asyncio.get_event_loop().run_in_executor(None, input, "")

    async def run(self):
        logger.info(f"🚀 Starting Ambeed Toxin Crawler (Sequential Warm-up Mode).")
        
        async with async_playwright() as p:
            browser, context = await self._init_browser(p)
            
            # 최초 접속
            page = await context.new_page()
            await page.goto("https://www.ambeed.com", timeout=60000)
            print("\n🛑 브라우저에서 차단(Cloudflare 등)이 있다면 풀고 엔터를 누르세요 (최초 1회)")
            await self.wait_for_user_input()
            await page.close()
            
            pages_in_batch = 0
            
            for page_num in range(self.page_start, self.max_pages_limit + 1):
                is_warmup = page_num < self.target_page
                # URL 패턴: adc-toxins.html?pagesize=20&pageindex=1
                target_url = f"{self.BASE_URL}?pagesize=20&pageindex={page_num}"
                
                logger.info(f"🎯 [Target] Processing Page {page_num}: {target_url}")
                
                # 배치 휴식 (5~6페이지쯤 막히는 것 대비 천천히)
                if pages_in_batch >= self.batch_size:
                    logger.info(f"🛑 Batch Resting ({self.rest_duration}s)... 천천히 진행합니다.")
                    await asyncio.sleep(self.rest_duration)
                    pages_in_batch = 0

                items_saved, items_found = await self.process_list_page(context, target_url, is_warmup=False)
                
                if items_found == 0:
                    logger.warning(f"⚠️ Page {page_num} blocked or empty. Please check browser.")
                    print("\n🛑 차단이 발생했거나 상품이 없습니다. 브라우저를 확인해주세요.")
                    print("계속하려면 엔터, 중단하려면 Ctrl+C를 누르세요.")
                    await self.wait_for_user_input("차단 해제 후 엔터...")
                    # 재시도
                    items_saved, items_found = await self.process_list_page(context, target_url, is_warmup=False)
                elif items_saved == 0:
                    logger.info(f"⏭️ Page {page_num}: All {items_found} items skipped (already in DB).")
                
                pages_in_batch += 1
                # 페이지 간 랜덤 딜레이
                await asyncio.sleep(random.uniform(10, 20))

                

                            if browser: await browser.close()

                            logger.info("✅ Crawler Finished.")

                

                    async def process_list_page(self, context, url, is_warmup):

                        page = await context.new_page()

                        items_processed = 0

                        total_links = 0

                        try:

                            await page.goto(url, wait_until="domcontentloaded", timeout=60000)

                            await asyncio.sleep(random.uniform(3, 6)) 

                            

                            product_links = await page.evaluate("""

                                () => {

                                    const anchors = Array.from(document.querySelectorAll('a'));

                                    return anchors

                                        .map(a => a.href)

                                        .filter(h => h.includes('/products/') && h.endsWith('.html'));

                                }

                            """)

                            

                            product_links = list(set(product_links))

                            total_links = len(product_links)

                            logger.info(f"🔎 Found {total_links} products on Page.")

                            

                            if not product_links:

                                return 0, 0

                

                            for p_url in product_links:

                                cat_no_match = re.search(r'/products/([a-zA-Z0-9-]+)\.html', p_url)

                                cat_no = cat_no_match.group(1) if cat_no_match else p_url.split('/')[-1].replace('.html', '')

                                

                                # 이미 DB에 있으면 스킵

                                if await self._exists_in_db(cat_no):

                                    # logger.info(f"Skipping {cat_no}")

                                    continue

                                

                                success = await self.process_product_detail(context, p_url)

                                if success:

                                    items_processed += 1

                                    # 상세 페이지 간 딜레이 (천천히)

                                    await asyncio.sleep(random.uniform(8, 15))

                

                        except Exception as e:

                            logger.error(f"❌ Page Error: {e}")

                            return 0, 0

                        finally:

                            await page.close()

                        

                        return items_processed, total_links

    async def _exists_in_db(self, cat_no: str) -> bool:
        try:
            res = supabase.table("commercial_reagents").select("id").eq("ambeed_cat_no", cat_no).limit(1).execute()
            return bool(res.data)
        except:
            return False

    async def process_product_detail(self, context, url):
        page = await context.new_page()
        try:
            await page.goto(url, wait_until="domcontentloaded", timeout=45000)
            
            data = await page.evaluate(
                """
                () => {
                    const res = { specs: {}, title: document.title };
                    const rows = document.querySelectorAll('table tr');
                    rows.forEach(row => {
                        const cells = row.querySelectorAll('td');
                        if (cells.length >= 2) {
                            let key = cells[0].innerText.trim().replace(/:$/, '');
                            let val = cells[1].innerText.trim();
                            
                            if (key.includes('CAS No')) key = 'cas_number';
                            else if (key.includes('Formula')) key = 'formula';
                            else if (key.includes('M.W')) key = 'molecular_weight';
                            else if (key.includes('SMILES Code')) key = 'smiles_code';
                            else if (key.includes('MDL No')) key = 'mdl_number';
                            else if (key.includes('InChI Key')) key = 'inchi_key';
                            else if (key.includes('Pubchem ID')) key = 'pubchem_id';
                            else if (key.includes('Synonyms')) key = 'synonyms';
                            
                            res.specs[key] = val;
                        }
                    });
                    return res;
                }
            """
            )
            
            specs = data['specs']
            cat_no_match = re.search(r'/products/([a-zA-Z0-9-]+)\.html', url)
            cat_no = cat_no_match.group(1) if cat_no_match else url.split('/')[-1].replace('.html', '')
            
            record = {
                "ambeed_cat_no": cat_no,
                "product_name": data['title'].split('|')[0].strip(),
                "category": "ADC Toxins", # 카테고리 설정
                "cas_number": specs.get('cas_number'),
                "formula": specs.get('formula'),
                "molecular_weight": specs.get('molecular_weight'),
                "smiles_code": specs.get('smiles_code'),
                "mdl_number": specs.get('mdl_number'),
                "source_name": "Ambeed",
                "product_url": url,
                "crawled_at": datetime.now().isoformat(),
                "properties": {
                    "inchi_key": specs.get('inchi_key'),
                    "pubchem_id": specs.get('pubchem_id'),
                    "synonyms": specs.get('synonyms') # Synonyms 저장
                }
            }
            supabase.table("commercial_reagents").upsert(record, on_conflict="ambeed_cat_no").execute()
            logger.info(f"✅ Saved: {cat_no} ({record['product_name']})")
            return True
        except Exception as e:
            logger.error(f"Failed to process {url}: {e}")
            return False
        finally:
            await page.close()

if __name__ == "__main__":
    crawler = AmbeedToxinSafeCrawler()
    asyncio.run(crawler.run())
