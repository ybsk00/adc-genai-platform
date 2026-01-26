import asyncio
from playwright.async_api import async_playwright
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DebugCBCategories")

async def debug():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        # Known Categories to verify
        urls = {
            "ADC Antibody": "https://www.creative-biolabs.com/adc/classify-adc-antibody-products-5.htm",
            "ADC Toxin": "https://www.creative-biolabs.com/adc/classify-adc-toxin-products-6.htm"
        }
        
        for cat, url in urls.items():
            logger.info(f"--- Checking {cat} ---")
            await page.goto(url)
            await asyncio.sleep(3)
            
            # Extract Pagination Link for Page 2
            page2_href = await page.evaluate("""
                () => {
                    const link = Array.from(document.querySelectorAll('.pages a, .pagination a')).find(a => a.innerText.trim() === '2');
                    return link ? link.href : null;
                }
            """)
            
            if page2_href:
                logger.info(f"✅ {cat} Page 2 URL: {page2_href}")
            else:
                logger.warning(f"⚠️ {cat} Page 2 link NOT found. Might be single page or different structure.")
                
        await browser.close()

if __name__ == "__main__":
    asyncio.run(debug())
