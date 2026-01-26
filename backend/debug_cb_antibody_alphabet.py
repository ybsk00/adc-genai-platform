import asyncio
from playwright.async_api import async_playwright
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DebugCBAlphabet")

async def debug():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        url = "https://www.creative-biolabs.com/adc/classify-adc-antibody-products-5.htm"
        logger.info(f"Visiting Antibody Page: {url}")
        await page.goto(url)
        await asyncio.sleep(5)
        
        # 1. Check for Alphabet Tabs (A-Z)
        logger.info("Scanning for Alphabet Tabs...")
        tabs = await page.evaluate("""
            () => {
                // Look for common alphabet tab containers
                const links = Array.from(document.querySelectorAll('.letter_list a, .zimu a, .alphabet a, ul.tabs a'));
                return links.map(a => ({text: a.innerText, href: a.href}));
            }
        """)
        
        if tabs:
            logger.info(f"✅ Found {len(tabs)} Alphabet Tabs: {[t['text'] for t in tabs[:5]]}...")
        else:
            logger.warning("⚠️ No explicit A-Z tabs found. Checking for 'Hot Targets' or Sub-menus.")
            
        # 2. Check content structure (Hot Targets / Combinations)
        content_headers = await page.evaluate("""
            () => Array.from(document.querySelectorAll('h2, h3, .title')).map(el => el.innerText)
        """)
        logger.info(f"Page Headers: {content_headers}")
        
        await browser.close()

if __name__ == "__main__":
    asyncio.run(debug())
