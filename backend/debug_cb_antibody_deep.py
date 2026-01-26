import asyncio
from playwright.async_api import async_playwright
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DebugCBAntibody")

async def debug():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        url = "https://www.creative-biolabs.com/adc/classify-adc-antibody-products-5.htm"
        logger.info(f"Visiting Antibody Page: {url}")
        await page.goto(url)
        await asyncio.sleep(5)
        
        # 1. Check for ANY pagination-like elements
        logger.info("Scanning for pagination controls...")
        pagination = await page.evaluate("""
            () => {
                const links = Array.from(document.querySelectorAll('.pages a, .pagination a, a.next, a.prev'));
                return links.map(a => ({text: a.innerText, href: a.href}));
            }
        """)
        
        if pagination:
            logger.info(f"Found Pagination Links: {pagination}")
        else:
            logger.warning("No standard pagination links found.")
            
        # 2. Check Item Count on Page 1
        count1 = await page.evaluate("document.querySelectorAll('.pro_list li, .list_pro li').length")
        item1 = await page.evaluate("""
            () => {
                const el = document.querySelector('.pro_list li a, .list_pro li a');
                return el ? el.innerText.trim() : "No Item";
            }
        """)
        logger.info(f"Page 1 Item Count: {count1}")
        logger.info(f"Page 1 First Item: {item1}")

        # 3. Force Visit Page 2
        url2 = "https://www.creative-biolabs.com/adc/classify-adc-antibodies-5.htm?page=2" # Try 'antibodies' pattern
        logger.info(f"Trying Constructed Page 2: {url2}")
        await page.goto(url2)
        await asyncio.sleep(5)
        
        item2 = await page.evaluate("""
            () => {
                const el = document.querySelector('.pro_list li a, .list_pro li a');
                return el ? el.innerText.trim() : "No Item";
            }
        """)
        
        if item1 == item2:
            logger.error("❌ Page 2 content is IDENTICAL to Page 1. Likely single page.")
        else:
            logger.info("✅ Page 2 content is DIFFERENT! Hidden pagination exists.")
            
        await browser.close()

if __name__ == "__main__":
    asyncio.run(debug())
