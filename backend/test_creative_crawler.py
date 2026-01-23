import asyncio
import logging
from app.services.creative_biolabs_crawler import creative_crawler

# Configure logging
logging.basicConfig(level=logging.INFO)

async def test_crawler():
    print("🚀 Testing Creative Biolabs Stealth Crawler...")
    
    # 1. Test PubChem API
    print("\n🧪 [Step 1] Testing PubChem API (Mylotarg CAS: 220572-57-0)...")
    smiles = await creative_crawler._fetch_pubchem_smiles("220572-57-0")
    print(f"   Result: {smiles}")
    if smiles:
        print("   ✅ PubChem API Passed")
    else:
        print("   ❌ PubChem API Failed")
        
    # 2. Test Gemini Enrichment
    print("\n🧪 [Step 2] Testing Gemini Enrichment...")
    desc = "Mylotarg is a CD33-directed antibody-drug conjugate (ADC)."
    data = await creative_crawler._enrich_with_gemini(desc)
    print(f"   Result: {data}")
    if data.get("target") == "CD33":
        print("   ✅ Gemini Enrichment Passed")
    else:
        print("   ❌ Gemini Enrichment Failed")

    # 3. Test Browser Launch (Stealth)
    print("\n🧪 [Step 3] Testing Playwright Stealth...")
    try:
        from playwright.async_api import async_playwright
        async with async_playwright() as p:
            context = await creative_crawler._init_browser(p)
            page = await context.new_page()
            await page.goto("https://bot.sannysoft.com/")
            await page.screenshot(path="stealth_check.png")
            print("   ✅ Browser Launched & Screenshot Saved (stealth_check.png)")
            await context.close()
    except Exception as e:
        print(f"   ❌ Browser Test Failed: {e}")

if __name__ == "__main__":
    asyncio.run(test_crawler())
