import asyncio
from playwright.async_api import async_playwright

async def debug():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        url = "https://www.creative-biolabs.com/adc/doxorubicin-716.htm"
        print(f"Navigating to {url}...")
        await page.goto(url)
        
        # Get Title
        title = await page.title()
        print(f"Title: {title}")
        
        # Get Content of interest
        content = await page.content()
        
        # Check for proshow-unit
        has_proshow = await page.evaluate("() => document.querySelector('.proshow-unit') !== null")
        print(f"Has .proshow-unit: {has_proshow}")
        
        # Dump a portion of HTML to file
        with open("backend/payload_debug.html", "w", encoding="utf-8") as f:
            f.write(content)
            
        await browser.close()

if __name__ == "__main__":
    asyncio.run(debug())
