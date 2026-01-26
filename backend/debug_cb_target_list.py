import asyncio
from playwright.async_api import async_playwright
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DebugCBTargetList")

async def debug():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        # Specific Target Page that failed
        url = "https://www.creative-biolabs.com/bsab/symbolsearch-ang2%2520%2526%2520vegfa.htm"
        logger.info(f"Visiting Target Page: {url}")
        await page.goto(url, timeout=60000)
        await asyncio.sleep(5)
        
        # Check for ANY links in content area
        links = await page.evaluate("""
            () => {
                // Try to find any link that looks like a product
                return Array.from(document.querySelectorAll('a'))
                    .filter(a => a.href.includes('/bsab/') && a.href.includes('-12')) // heuristic for product IDs
                    .map(a => ({text: a.innerText, href: a.href}))
                    .slice(0, 10);
            }
        """)
        
        logger.info(f"Potential Product Links: {links}")
        
        # Check DOM structure for lists
        structure = await page.evaluate("""
            () => {
                const lists = document.querySelectorAll('ul, ol, div.list');
                return Array.from(lists).map(l => l.className);
            }
        """)
        logger.info(f"List Classes Found: {structure}")
        
        await browser.close()

if __name__ == "__main__":
    asyncio.run(debug())
