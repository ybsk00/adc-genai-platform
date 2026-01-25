import asyncio
from playwright.async_api import async_playwright
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("TestRealPagination")

async def test():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        # Test Page 1
        url1 = "https://www.ambeed.com/adc-toxins.html?pagesize=20&pageindex=1"
        logger.info(f"Visiting P1: {url1}")
        await page.goto(url1)
        await asyncio.sleep(3)
        p1_items = await page.evaluate("() => Array.from(document.querySelectorAll('a[href*=\"/products/\"]')).map(a => a.href).slice(0, 5)")
        logger.info(f"P1 Sample Items: {p1_items}")
        
        # Test Page 2
        url2 = "https://www.ambeed.com/adc-toxins.html?pagesize=20&pageindex=2"
        logger.info(f"Visiting P2: {url2}")
        await page.goto(url2)
        await asyncio.sleep(3)
        p2_items = await page.evaluate("() => Array.from(document.querySelectorAll('a[href*=\"/products/\"]')).map(a => a.href).slice(0, 5)")
        logger.info(f"P2 Sample Items: {p2_items}")
        
        if p1_items == p2_items:
            logger.error("❌ CONTENT IS THE SAME! pageindex=2 IS IGNORED.")
            
            # Try alternative: page=2
            url_alt = "https://www.ambeed.com/adc-toxins.html?page=2"
            logger.info(f"Trying Alternative: {url_alt}")
            await page.goto(url_alt)
            await asyncio.sleep(3)
            p_alt_items = await page.evaluate("() => Array.from(document.querySelectorAll('a[href*=\"/products/\"]')).map(a => a.href).slice(0, 5)")
            logger.info(f"Alt Sample Items: {p_alt_items}")
            if p1_items != p_alt_items:
                logger.info("✅ Alternative page=2 WORKS!")
        else:
            logger.info("✅ pageindex=2 WORKS! Content is different.")
            
        await browser.close()

if __name__ == "__main__":
    asyncio.run(test())
