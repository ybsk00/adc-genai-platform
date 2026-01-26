import asyncio
from playwright.async_api import async_playwright
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DebugCBPagination")

async def debug():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        # Creative Biolabs ADC Linker Category
        base_url = "https://www.creative-biolabs.com/adc/classify-adc-linker-products-7.htm"
        logger.info(f"Visiting {base_url}...")
        await page.goto(base_url)
        await asyncio.sleep(5)
        
        # 1. Extract Pagination Links (Onclick)
        logger.info("Extracting pagination info (onclick)...")
        pagination = await page.evaluate("""
            () => {
                const pageLinks = Array.from(document.querySelectorAll('.pages a, .pagination a')).map(a => ({
                    text: a.innerText, 
                    href: a.href,
                    onclick: a.getAttribute('onclick')
                }));
                return pageLinks;
            }
        """)
        
        for link in pagination[:5]:
            logger.info(f"Link: {link}")

        # 2. Check Item 1 on Page 1 (Try broader selectors)
        item1 = await page.evaluate("""
            () => {
                // Try multiple selectors for product list
                const el = document.querySelector('.pro_list li a, .list_pro li a, .prod_list li a, .cp_list li a');
                return el ? el.innerText.trim() : "No Item Found";
            }
        """)
        logger.info(f"Page 1 First Item: {item1[:50]}...")
        
        # 3. Test ?page=2
        target_page_2 = f"{base_url}?page=2"
        logger.info(f"Visiting constructed Page 2: {target_page_2}")
        await page.goto(target_page_2)
        await asyncio.sleep(5)
        
        item2 = await page.evaluate("""
            () => {
                const el = document.querySelector('.pro_list .pro_name, .list_pro .name, .prod_list h3 a');
                return el ? el.innerText : "No Item Found";
            }
        """)
        logger.info(f"Page 2 First Item: {item2[:50]}...")
        
        if item1 == item2:
            logger.error("⚠️ Page 1 and Page 2 content are IDENTICAL! ?page=2 might be wrong.")
        else:
            logger.info("✅ Page 2 is different. ?page=2 works.")
            
        await browser.close()

if __name__ == "__main__":
    asyncio.run(debug())
