import asyncio
from playwright.async_api import async_playwright

async def debug_html():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        url = "https://www.ambeed.com/search?keyword=ADC+Toxin"
        print(f"Navigating to {url}...")
        try:
            await page.goto(url, wait_until="domcontentloaded", timeout=60000)
            content = await page.content()
            
            # Save to file
            with open("ambeed_debug.html", "w", encoding="utf-8") as f:
                f.write(content)
            print("Saved HTML to ambeed_debug.html")
            
            # Print some links to see what's there
            links = await page.evaluate("""
                Array.from(document.querySelectorAll('a'))
                    .map(a => a.href)
                    .filter(href => href.includes('product'))
            """)
            print("Found product-like links:", links[:5])
            
        except Exception as e:
            print(f"Error: {e}")
        finally:
            await browser.close()

if __name__ == "__main__":
    asyncio.run(debug_html())
