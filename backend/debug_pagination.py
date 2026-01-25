import asyncio
from playwright.async_api import async_playwright
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DebugPagination")

async def debug():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        url = "https://www.ambeed.com/adc-toxins.html"
        logger.info(f"Visiting {url}...")
        await page.goto(url)
        await asyncio.sleep(3)
        
        # 1. Extract Pagination Links
        logger.info("Extracting pagination info...")
        pagination = await page.evaluate("""
            () => {
                const nextBtn = document.querySelector('a[rel="next"], a.next, .pages a.next');
                const pageLinks = Array.from(document.querySelectorAll('.pages a, .pagination a')).map(a => ({text: a.innerText, href: a.href}));
                return {
                    next_href: nextBtn ? nextBtn.href : null,
                    all_links: pageLinks
                };
            }
        """)
        
        logger.info(f"Next Button HREF: {pagination['next_href']}")
        logger.info(f"All Page Links: {pagination['all_links']}")
        
        # 2. Check Item 1 on Page 1
        item1 = await page.evaluate("document.querySelector('.product-item, .item').innerText")
        logger.info(f"Page 1 First Item: {item1[:50]}...")
        
        # 3. Visit Page 2 (Constructed or Extracted)
        target_page_2 = "https://www.ambeed.com/adc-toxins.html?page=2"
        logger.info(f"Visiting constructed Page 2: {target_page_2}")
        await page.goto(target_page_2)
        await asyncio.sleep(3)
        
        item2 = await page.evaluate("document.querySelector('.product-item, .item').innerText")
        logger.info(f"Page 2 First Item: {item2[:50]}...")
        
        if item1 == item2:
            logger.error("⚠️ Page 1 and Page 2 content are IDENTICAL! URL param ?page=2 is ignored.")
        else:
            logger.info("✅ Page 2 is different.")
            
        await browser.close()

if __name__ == "__main__":
    asyncio.run(debug())
