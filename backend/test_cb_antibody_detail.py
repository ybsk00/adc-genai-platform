import asyncio
from playwright.async_api import async_playwright
import logging
import json

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("TestCBAntibodyDetail")

async def test():
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        
        # 1. Visit By Targets Page
        base_url = "https://www.creative-biolabs.com/bsab/category/by-targets-1377.htm"
        logger.info(f"Step 1: Visiting {base_url}")
        await page.goto(base_url, timeout=60000)
        
        # 2. Simulate Alphabet Click (e.g., 'A') - Testing logical navigation
        # In reality, we might just grab all symbolsearch links if they are all loaded or loadable.
        # But for this test, we'll try to find a known target link directly.
        
        target_href = "symbolsearch-ang2%2520%2526%2520vegfa.htm" # Ang2 & VEGFA
        # Fix URL encoding for finding the link
        target_selector = f'a[href*="{target_href}"]'
        
        logger.info(f"Step 2: Looking for Target Link: {target_href}")
        # Note: In real crawl, we iterate alphabet. Here we assume the link might be visible or we go directly.
        # Let's try direct navigation to target page to simulate "found link".
        target_full_url = f"https://www.creative-biolabs.com/bsab/{target_href}"
        await page.goto(target_full_url, timeout=60000)
        await asyncio.sleep(3)
        
        # 3. Find Product Link
        logger.info("Step 3: Finding Product Link...")
        product_link = await page.evaluate("""
            () => {
                const link = document.querySelector('.pro_list li a, .list_pro li a');
                return link ? link.href : null;
            }
        """)
        
        if not product_link:
            logger.error("❌ No product link found on target page.")
            await browser.close()
            return

        logger.info(f"✅ Found Product: {product_link}")
        await page.goto(product_link, timeout=60000)
        await asyncio.sleep(3)
        
        # 4. Extract Detail Data
        logger.info("Step 4: Extracting Specifications...")
        data = await page.evaluate(r"""
            () => {
                const specs = {};
                // Extract Tables
                document.querySelectorAll('table tr').forEach(row => {
                    const cells = row.querySelectorAll('th, td');
                    if (cells.length >= 2) {
                        const key = cells[0].innerText.trim().replace(/:$/, '');
                        const val = cells[1].innerText.trim();
                        specs[key] = val;
                    }
                });
                
                // Extract "Targets" section specifically if structured differently
                const targetInfo = {};
                // Custom logic for Gene ID / UniProt ID which might be in a separate text block or table
                // Looking for text patterns in the whole body or specific sections
                const bodyText = document.body.innerText;
                const geneIdMatch = bodyText.match(/Gene ID\s*(\d+)/i);
                const uniprotMatch = bodyText.match(/UniProt ID\s*([A-Z0-9]+)/i);
                
                if (geneIdMatch) targetInfo['Gene ID'] = geneIdMatch[1];
                if (uniprotMatch) targetInfo['UniProt ID'] = uniprotMatch[1];
                
                return { specs, targetInfo };
            }
        """)
        
        logger.info(f"✅ Extracted Data:\n{json.dumps(data, indent=2)}")
        
        await browser.close()

if __name__ == "__main__":
    asyncio.run(test())
