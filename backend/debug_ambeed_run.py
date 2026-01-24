import asyncio
from playwright.async_api import async_playwright
from fake_useragent import UserAgent

async def debug_ambeed():
    ua = UserAgent()
    user_agent = ua.random
    print(f"Using User-Agent: {user_agent}")

    urls = [
        "https://www.ambeed.com/adc-toxins.html",
        "https://www.ambeed.com/antibody-drug-conjugate-adc-related.html",
        "https://www.ambeed.com/payload-linker-conjugate.html"
    ]
    
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        context = await browser.new_context(
            user_agent=user_agent,
            viewport={'width': 1920, 'height': 1080}
        )
        page = await context.new_page()
        
        # urls already defined above
        
        for url in urls:
            print(f"\n--- Testing Category: {url} ---")
            try:
                await page.goto(url, wait_until="domcontentloaded", timeout=60000)
                await page.wait_for_timeout(5000)
                
                # Extract product links
                links = await page.evaluate("""
                    Array.from(document.querySelectorAll('a[href*="/products/"], a[href*="/record/"], .product-item a, .pro-list a'))
                        .map(a => a.href)
                        .filter(href => href.includes('/products/') || href.includes('/record/'))
                """)
                links = list(set(links))
                print(f"Found {len(links)} product links.")
                
                if links:
                    product_url = links[0]
                    print(f"Testing Product Page: {product_url}")
                    await page.goto(product_url, wait_until="domcontentloaded", timeout=60000)
                    await page.wait_for_timeout(3000)
                    
                    # Screenshot
                    await page.screenshot(path="ambeed_product_debug.png")
                    print("Saved product screenshot to ambeed_product_debug.png")
                    
                    # Save HTML
                    content = await page.content()
                    with open("ambeed_product_debug.html", "w", encoding="utf-8") as f:
                        f.write(content)
                    print("Saved product HTML to ambeed_product_debug.html")

                    # Test Extraction Logic
                    title = await page.title()
                    print(f"Product Title: {title}")
                    
                    # Dump Body Text
                    body_text = await page.inner_text("body")
                    with open("ambeed_body_text.txt", "w", encoding="utf-8") as f:
                        f.write(body_text)
                    print("Saved body text to ambeed_body_text.txt")

                    # JSON-LD Extraction
                    json_ld = await page.evaluate("""
                        () => {
                            const scripts = Array.from(document.querySelectorAll('script[type="application/ld+json"]'));
                            return scripts.map(s => JSON.parse(s.innerText));
                        }
                    """)
                    print(f"JSON-LD Data: {json_ld}")

                    async def extract_by_label(label):
                        return await page.evaluate(f"""
                            (label) => {{
                                // Improved extraction strategy
                                const elements = Array.from(document.querySelectorAll('td, th, div, span, p, b, strong, li'));
                                
                                for (const el of elements) {{
                                    // STRICT CHECK: Skip elements that contain many children or too much text
                                    // We want leaf nodes or nodes close to leaves
                                    if (el.children.length > 2 || el.innerText.length > 300) continue;
                                    
                                    const text = el.innerText.trim();
                                    if (text.includes(label)) {{
                                        
                                        // Strategy 1: "Label: Value" in same element
                                        // Use Regex to be safer
                                        const safeLabel = label.replace(/[.*+?^${{}}()|[\\]\\\\]/g, '\\\\$&');
                                        const regex = new RegExp(safeLabel + '\\\\s*[:\\\\uff1a]\\\\s*([^\\\\n\\\\r]+)', 'i');
                                        const match = text.match(regex);
                                        
                                        if (match) {{
                                            return match[1].trim();
                                        }}
                                        
                                        // Strategy 2: Value in next sibling (Label is in <b> or <td>, value in next node)
                                        // Strict check: element text should be mostly just the label
                                        if (text.replace(':', '').trim() === label || text.length < label.length + 5) {{
                                            if (el.nextElementSibling) {{
                                                return el.nextElementSibling.innerText.trim();
                                            }}
                                            if (el.parentElement && el.parentElement.nextElementSibling) {{
                                                // Check if parent is small enough to be just a wrapper
                                                if (el.parentElement.innerText.length < 100) {{
                                                     return el.parentElement.nextElementSibling.innerText.trim();
                                                }}
                                            }}
                                        }}
                                    }}
                                }}
                                return null;
                            }}
                        """, label)

                    cat_no = await extract_by_label("Catalog No")
                    cas_no = await extract_by_label("CAS No")
                    formula = await extract_by_label("Molecular Formula")
                    if not formula: formula = await extract_by_label("Formula")
                    mw = await extract_by_label("Molecular Weight")
                    if not mw: mw = await extract_by_label("M.W")
                    
                    smiles = await extract_by_label("SMILES Code")
                    if not smiles:
                        smiles = await extract_by_label("SMILES")

                    # Text Pattern Fallback
                    if not smiles:
                        smiles_from_text = await page.evaluate("""
                            () => {
                                const bodyText = document.body.innerText;
                                const match = bodyText.match(/SMILES Code\\s*:\\s*([^\\n]+)/);
                                if (match) return match[1].trim();
                                return null;
                            }
                        """)
                        if smiles_from_text:
                            smiles = smiles_from_text
                            print(f"Found SMILES via Text Pattern: {smiles}")

                    print(f"Extracted Data:")
                    print(f" - Catalog No: {cat_no}")
                    print(f" - CAS No: {cas_no}")
                    print(f" - Formula: {formula}")
                    print(f" - MW: {mw}")
                    print(f" - SMILES: {smiles}")
                    
            except Exception as e:
                print(f"Error testing {url}: {e}")
        
        await browser.close()

if __name__ == "__main__":
    asyncio.run(debug_ambeed())
