import asyncio
from playwright.async_api import async_playwright
import logging
from urllib.parse import urljoin

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("CB_Structure_Analysis")

BASE_URL = "https://www.creative-biolabs.com/bsab/"

async def analyze():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        # 1. Analyze 'By Targets' Page
        target_url = "https://www.creative-biolabs.com/bsab/category/by-targets-1377.htm"
        logger.info(f"🔍 Analyzing Targets Page: {target_url}")
        await page.goto(target_url, timeout=60000)
        await asyncio.sleep(5)
        
        target_links = await page.evaluate("""
            () => {
                return Array.from(document.querySelectorAll('a[href*="symbolsearch-"]'))
                    .map(a => a.href)
            }
        """)
        unique_targets = list(set(target_links))
        logger.info(f"✅ Found {len(unique_targets)} Unique Targets (e.g., {unique_targets[:3]})")
        
        # 2. Analyze a Sample Target Page (Deep Dive)
        if unique_targets:
            sample_target = unique_targets[0]
            logger.info(f"🔍 Deep Diving into Sample Target: {sample_target}")
            await page.goto(sample_target, timeout=60000)
            await asyncio.sleep(5)
            
            # Check for Products
            products = await page.evaluate("""
                () => document.querySelectorAll('.pro_list li, .list_pro li').length
            """)
            logger.info(f"📦 Sample Target contains {products} products.")
            
            # Check Pagination on Target Page
            pagination = await page.evaluate("""
                () => document.querySelector('.pages, .pagination') ? true : false
            """)
            logger.info(f"📄 Pagination present on Target Page? {pagination}")

        # 3. Analyze Main Categories (Linker, Toxin, etc. from ADC section)
        # Assuming ADC section URLs from previous knowledge
        adc_base = "https://www.creative-biolabs.com/adc/"
        logger.info(f"🔍 Checking ADC Main Page: {adc_base}")
        await page.goto(adc_base, timeout=60000)
        
        adc_categories = await page.evaluate("""
            () => Array.from(document.querySelectorAll('a[href*="classify-"]'))
                .map(a => ({text: a.innerText, href: a.href}))
        """)
        logger.info(f"✅ Found {len(adc_categories)} ADC Categories: {[c['text'] for c in adc_categories[:5]]}")

        await browser.close()

if __name__ == "__main__":
    asyncio.run(analyze())
