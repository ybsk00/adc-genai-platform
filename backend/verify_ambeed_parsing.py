import asyncio
import logging
from playwright.async_api import async_playwright
from fake_useragent import UserAgent

# Logging setup
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

async def verify_ambeed_parsing():
    """
    Verify SMILES extraction from Ambeed product page without using AI.
    Target URL: A known Ambeed product page with SMILES.
    """
    # Example URL (MMAE) - Replace with a valid one if needed
    # Using a search result or a specific product page if known. 
    # Since I don't have a specific URL guaranteed to be alive, I'll use a search and click first result or use a hardcoded one if I knew it.
    # Let's try a direct product page if possible, or search for 'MMAE'.
    
    target_url = "https://www.ambeed.com/search?keyword=MMAE" 
    
    ua = UserAgent()
    user_agent = ua.random
    
    logger.info(f"🚀 Starting Ambeed Parsing Verification...")
    logger.info(f"Target URL: {target_url}")

    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        context = await browser.new_context(
            user_agent=user_agent,
            viewport={'width': 1920, 'height': 1080}
        )
        
        # Stealth
        await context.add_init_script("""
            Object.defineProperty(navigator, 'webdriver', {
                get: () => undefined
            });
        """)
        
        page = await context.new_page()
        
        try:
            # 1. Go to Search Page
            await page.goto(target_url, wait_until="networkidle", timeout=60000)
            logger.info("✅ Search page loaded.")
            
            # 2. Click first product
            # Selector from ambeed_crawler.py
            product_link = await page.evaluate("""
                () => {
                    const a = document.querySelector('a[href*="/products/"], a[href*="/record/"]');
                    return a ? a.href : null;
                }
            """)
            
            if not product_link:
                logger.error("❌ No product link found on search page.")
                return
                
            logger.info(f"🔗 Found product link: {product_link}")
            await page.goto(product_link, wait_until="domcontentloaded")
            logger.info("✅ Product page loaded.")
            
            # 3. Extract Data (Logic from ambeed_crawler.py)
            
            # Helper to extract by label
            async def extract_by_label(label):
                return await page.evaluate(f"""
                    (label) => {{
                        const elements = Array.from(document.querySelectorAll('td, th, div, span, p'));
                        for (const el of elements) {{
                            if (el.innerText.includes(label)) {{
                                const text = el.innerText;
                                if (text.includes(':')) return text.split(':')[1].trim();
                                if (el.nextElementSibling) return el.nextElementSibling.innerText.trim();
                            }}
                        }}
                        return null;
                    }}
                """, label)

            cat_no = await extract_by_label("Catalog No")
            cas_no = await extract_by_label("CAS No")
            smiles = await extract_by_label("SMILES Code")
            
            if not smiles:
                 smiles = await extract_by_label("SMILES")
            
            # Text pattern fallback
            if not smiles:
                smiles_from_text = await page.evaluate("""
                    () => {
                        const bodyText = document.body.innerText;
                        const match = bodyText.match(/SMILES Code\\s*:\\s*([A-Za-z0-9@#\\(\\)\\[\\]\\/\\\\=+-]+)/);
                        return match ? match[1].trim() : null;
                    }
                """)
                if smiles_from_text:
                    smiles = smiles_from_text
                    logger.info("✅ Found SMILES via Text Pattern.")

            logger.info("-" * 30)
            logger.info(f"📌 Extraction Results:")
            logger.info(f"   Catalog No: {cat_no}")
            logger.info(f"   CAS No:     {cas_no}")
            logger.info(f"   SMILES:     {smiles}")
            logger.info("-" * 30)
            
            if smiles:
                logger.info("✨ SUCCESS: SMILES extracted without AI!")
            else:
                logger.warning("⚠️ FAILURE: SMILES not found.")

        except Exception as e:
            logger.error(f"❌ Error: {e}")
        finally:
            await browser.close()

if __name__ == "__main__":
    asyncio.run(verify_ambeed_parsing())
