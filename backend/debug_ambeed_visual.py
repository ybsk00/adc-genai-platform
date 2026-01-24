import asyncio
import logging
from playwright.async_api import async_playwright

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("AMBEED_DEBUG")

async def debug_visual_check():
    url = "https://www.ambeed.com/adc-toxins.html"
    print(f"🚀 Navigating to: {url}")

    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True) # Change to False if you want to see it locally (requires GUI)
        context = await browser.new_context(
            user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
        )
        page = await context.new_page()

        try:
            await page.goto(url, wait_until="networkidle", timeout=60000)
            
            title = await page.title()
            print(f"✅ Page Title: {title}")
            
            # Screenshot for evidence
            await page.screenshot(path="ambeed_debug_snapshot.png", full_page=True)
            print("📸 Screenshot saved to: ambeed_debug_snapshot.png")
            
            # Count products
            count = await page.evaluate("""
                Array.from(document.querySelectorAll('a[href*="/products/"], a[href*="/record/"], .product-item a, .pro-list a'))
                    .map(a => a.href)
                    .filter(href => href.includes('/products/') || href.includes('/record/'))
                    .length
            """)
            print(f"🔢 Product Links Found: {count}")
            
            if count == 0:
                print("⚠️ No products found! Selectors might need update.")
                print("   Dumping HTML for inspection...")
                content = await page.content()
                with open("ambeed_debug_dump.html", "w", encoding="utf-8") as f:
                    f.write(content)
            else:
                print("✅ Parsing Logic seems OK.")

        except Exception as e:
            print(f"❌ Error during debug: {e}")
        finally:
            await browser.close()

if __name__ == "__main__":
    asyncio.run(debug_visual_check())
