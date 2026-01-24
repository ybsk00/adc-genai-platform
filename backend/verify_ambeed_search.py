import asyncio
from playwright.async_api import async_playwright

async def verify_url():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        url = "https://www.ambeed.com/search/Search.html?keyword=ADC+Toxin"
        print(f"Testing {url}...")
        try:
            response = await page.goto(url, wait_until="domcontentloaded", timeout=30000)
            print(f"  Status: {response.status}")
            title = await page.title()
            print(f"  Title: {title}")
            
            # Check for products
            # Ambeed search results usually have product links
            links = await page.evaluate("""
                Array.from(document.querySelectorAll('a'))
                    .map(a => a.href)
                    .filter(href => href.includes('/products/'))
            """)
            print(f"  Found {len(links)} product links.")
            if links:
                print(f"  Example link: {links[0]}")
                
        except Exception as e:
            print(f"  Error: {e}")
        
        await browser.close()

if __name__ == "__main__":
    asyncio.run(verify_url())
