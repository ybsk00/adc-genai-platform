import asyncio
import logging
import os
import json
from playwright.async_api import async_playwright
from dotenv import load_dotenv
from supabase import create_client

# Load Env
load_dotenv()

# Setup Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("TestOne")

# Init Supabase
url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_KEY")
if not url or not key:
    logger.error("Supabase credentials missing!")
    exit(1)

supabase = create_client(url, key)

async def test_one_product_db():
    # Target URL
    url = "https://www.ambeed.com/products/1807534-78-6.html" 
    
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        context = await browser.new_context()
        page = await context.new_page()
        
        logger.info(f"Visiting {url}...")
        await page.goto(url, timeout=60000)
        
        extracted = await page.evaluate("""
            () => {
                const data = {
                    cas_no: null,
                    formula: null,
                    mw: null,
                    smiles: null,
                    mdl_no: null,
                    product_name: document.title.split('|')[0].trim()
                };
                
                const clean = (text) => text ? text.replace(/:/g, '').trim() : null;

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
                return data;
            }
        """)
        
        # Prepare DB Data
        cat_no = "1807534-78-6" 
        
        # Move mdl_number to properties to avoid Schema Error
        db_item = {
            "ambeed_cat_no": cat_no,
            "product_name": extracted['product_name'],
            "product_url": url,
            "category": "Test Category",
            "cas_number": extracted['cas_no'],
            "formula": extracted['formula'],
            "molecular_weight": extracted['mw'],
            "smiles_code": extracted['smiles'],
            "source_name": "Ambeed",
            "crawled_at": "2026-01-26T00:00:00Z",
            "properties": {
                "mdl_number": extracted['mdl_no']
            }
        }
        
        logger.info(f"Attempting to save to DB (mdl_number inside properties): {db_item['ambeed_cat_no']}")
        
        try:
            res = supabase.table("commercial_reagents").upsert([db_item], on_conflict="ambeed_cat_no").execute()
            logger.info(f"DB Response: {res}")
            
            if res.data:
                logger.info("✅ Save SUCCESS!")
            else:
                logger.error("❌ Save FAILED (No data returned).")
                
        except Exception as e:
            logger.error(f"❌ DB Error: {e}")

        await browser.close()

if __name__ == "__main__":
    asyncio.run(test_one_product_db())