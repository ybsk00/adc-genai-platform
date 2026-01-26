import asyncio
from playwright.async_api import async_playwright

async def dump():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        await page.goto("https://www.creative-biolabs.com/bsab/category/by-targets-1377.htm")
        await asyncio.sleep(5)
        
        content = await page.content()
        with open("cb_bsab_dump.html", "w", encoding="utf-8") as f:
            f.write(content)
        print("✅ HTML saved to cb_bsab_dump.html")
        await browser.close()

if __name__ == "__main__":
    asyncio.run(dump())
