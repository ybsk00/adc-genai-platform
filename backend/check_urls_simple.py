import asyncio
from playwright.async_api import async_playwright

async def test_urls():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        urls = [
            "https://www.ambeed.com/search.html?keyword=ADC+Toxin",
            "https://www.ambeed.com/search-products.html?keyword=ADC+Toxin",
            "https://www.ambeed.com/search/keyword.html?keyword=ADC+Toxin",
            "https://www.ambeed.com/search-products?keyword=ADC+Toxin",
            "https://www.ambeed.com/products/search?keyword=ADC+Toxin"
        ]
        
        for url in urls:
            try:
                response = await page.goto(url, wait_until="domcontentloaded", timeout=10000)
                title = await page.title()
                if "404" not in title and response.status == 200:
                    print(f"SUCCESS: {url}")
                    return
            except:
                pass
        
        print("NONE_FOUND")
        await browser.close()

if __name__ == "__main__":
    asyncio.run(test_urls())
