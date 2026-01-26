import asyncio
from playwright.async_api import async_playwright
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DebugBsAb")

async def debug():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        # New URL provided by user
        url = "https://www.creative-biolabs.com/bsab/category/by-targets-1377.htm"
        logger.info(f"Visiting BsAb Target Page: {url}")
        await page.goto(url)
        await asyncio.sleep(5)
        
        # 1. Check for Alphabet Tabs (A-Z)
        logger.info("Scanning for Alphabet Tabs...")
        tabs = await page.evaluate("""
            () => {
                // Common class for alphabet filters
                const links = Array.from(document.querySelectorAll('.character a, .zimu a, .alphabet a'));
                return links.map(a => ({text: a.innerText, href: a.href}));
            }
        """)
        
        if tabs:
            logger.info(f"✅ Found {len(tabs)} Alphabet Tabs: {[t['text'] for t in tabs[:5]]}...")
        else:
            logger.warning("⚠️ Alphabet tabs not found with standard selectors.")
            
        # 2. Check for Target Links
        targets = await page.evaluate("""
            () => {
                const links = Array.from(document.querySelectorAll('.list_pro li a, .pro_list li a'));
                return links.map(a => ({text: a.innerText, href: a.href})).slice(0, 5);
            }
        """)
        logger.info(f"Sample Targets: {targets}")
        
        await browser.close()

if __name__ == "__main__":
    asyncio.run(debug())
