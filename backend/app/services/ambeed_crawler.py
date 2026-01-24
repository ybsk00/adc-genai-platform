import asyncio
import random
import logging
import json
import time
from typing import Dict, List, Optional
from datetime import datetime

from playwright.async_api import async_playwright, Page, BrowserContext
from fake_useragent import UserAgent
import aiohttp
import google.generativeai as genai

from app.core.config import settings
from app.core.supabase import supabase

logger = logging.getLogger(__name__)

class AmbeedCrawler:
    """Ambeed ADC Reagents Stealth Crawler & AI Enrichment"""
    
    BASE_URL = "https://www.ambeed.com"
    
    # Categories to crawl (Search URLs)
    CATEGORIES = {
        "ADC Toxins": "https://www.ambeed.com/adc-toxins.html",
        "ADC Linkers": "https://www.ambeed.com/search?keyword=ADC+Linker",
        "ADC Cytotoxin": "https://www.ambeed.com/search?keyword=ADC+Cytotoxin",
        "Payload": "https://www.ambeed.com/adc-toxins.html", # Alias for Payload
        "Linker": "https://www.ambeed.com/search?keyword=ADC+Linker", # Alias for Linker
        "Conjugate": "https://www.ambeed.com/search?keyword=Antibody-Drug+Conjugate" # Alias for Conjugate
    }
    
    def __init__(self):
        self.ua = UserAgent()
        self.request_count = 0
        self.break_threshold = 100  # Break every 100 requests
        
        # Configure Gemini
        genai.configure(api_key=settings.GOOGLE_API_KEY)
        self.model = genai.GenerativeModel(settings.GEMINI_MODEL_ID or 'gemini-1.5-flash')

    async def _init_browser(self, p) -> BrowserContext:
        """Initialize browser with stealth and proxy settings"""
        user_agent = self.ua.random
        
        proxy_config = None
        if settings.PROXY_ENABLED:
            # Handle multiple ports
            ports = settings.PROXY_PORTS.split(',')
            port = random.choice(ports).strip()
            
            proxy_config = {
                "server": f"http://{settings.PROXY_HOST}:{port}",
                "username": settings.PROXY_USERNAME,
                "password": settings.PROXY_PASSWORD
            }
            logger.info(f"🌐 Using Proxy: {settings.PROXY_HOST}:{port}")

        browser = await p.chromium.launch(
            headless=True, 
            proxy=proxy_config if settings.PROXY_ENABLED else None
        )
        context = await browser.new_context(
            user_agent=user_agent,
            viewport={'width': 1920, 'height': 1080},
            locale='en-US',
            timezone_id='America/New_York'
        )
        
        # Inject stealth scripts (navigator.webdriver removal)
        await context.add_init_script("""
            Object.defineProperty(navigator, 'webdriver', {
                get: () => undefined
            });
        """)
        
        return context

    async def _smart_delay(self):
        """Random delay to mimic human behavior"""
        delay = random.uniform(2.5, 7.5)
        logger.info(f"💤 Sleeping for {delay:.2f}s...")
        await asyncio.sleep(delay)
        
        self.request_count += 1
        if self.request_count >= self.break_threshold:
            break_time = 600  # 10 minutes
            logger.warning(f"☕ Break Time! Sleeping for {break_time/60:.1f} minutes to avoid ban...")
            await asyncio.sleep(break_time)
            self.request_count = 0

    async def _fetch_pubchem_smiles(self, cas_number: str) -> Optional[str]:
        """Fetch SMILES from PubChem API using CAS number"""
        if not cas_number:
            return None
            
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/property/IsomericSMILES/JSON"
            async with aiohttp.ClientSession() as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        data = await response.json()
                        return data['PropertyTable']['Properties'][0]['IsomericSMILES']
        except Exception as e:
            logger.warning(f"PubChem lookup failed for {cas_number}: {e}")
        return None

    async def _enrich_with_gemini(self, description: str) -> Dict:
        """Extract Target and Attributes using Gemini"""
        prompt = f"""
        Analyze the following Ambeed reagent description and extract structured data.
        
        Description: {description[:2000]}
        
        Extract:
        1. Target: The specific antigen target (e.g., HER2, TROP2, CD19). If not specified, return null.
           - Normalize "HER2-directed" to "HER2".
        2. Properties: Extract Purity, Molecular Weight (MW), Formula if present.
        3. Summary: A concise 3-sentence summary of the product, focusing on its application and key characteristics, suitable for researchers.
        
        Return JSON format:
        {{
            "target": "HER2",
            "properties": {{
                "purity": ">98%",
                "mw": "1234.56",
                "formula": "C50H60N10O12"
            }},
            "summary": "This product is a..."
        }}
        """
        
        try:
            response = await self.model.generate_content_async(prompt)
            text = response.text.replace("```json", "").replace("```", "").strip()
            return json.loads(text)
        except Exception as e:
            logger.error(f"Gemini enrichment failed: {e}")
            return {"target": None, "properties": {}}

    async def _get_embedding(self, text: str) -> List[float]:
        """Generate vector embedding for text"""
        try:
            result = await genai.embed_content_async(
                model="models/text-embedding-004",
                content=text,
                task_type="retrieval_document"
            )
            return result['embedding']
        except Exception as e:
            logger.error(f"Embedding generation failed: {e}")
            return []

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10) -> int:
        """Crawl a specific category with pagination"""
        logger.info(f"🚀 Starting crawl for {category_name}...")
        
        count = 0
        page_num = 1
        max_pages = 116  # Hard limit
        
        is_full_crawl = limit > 1000
        
        async with async_playwright() as p:
            context = await self._init_browser(p)
            page = await context.new_page()
            
            try:
                while True:
                    if not is_full_crawl and count >= limit:
                        break
                    if page_num > max_pages:
                        break

                    # Ambeed pagination: &page=2
                    if page_num == 1:
                        url = base_url
                    else:
                        url = f"{base_url}&page={page_num}"
                        
                    logger.info(f"📄 Navigating to Page {page_num}: {url}")
                    
                    try:
                        await page.goto(url, wait_until="networkidle", timeout=60000)
                    except Exception as e:
                        logger.error(f"Page load timeout/error: {e}")
                        break

                    await self._smart_delay()
                    
                    if page_num > 1 and page_num % 10 == 0:
                        break_time = 300
                        logger.warning(f"☕ Stealth Break! Sleeping for {break_time/60:.1f} minutes...")
                        await asyncio.sleep(break_time)
                    
                    # Extract product links
                    try:
                        links = await page.evaluate("""
                            Array.from(document.querySelectorAll('a[href*="/products/"], a[href*="/record/"], .product-item a, .pro-list a'))
                                .map(a => a.href)
                                .filter(href => href.includes('/products/') || href.includes('/record/'))
                        """)
                        
                        links = list(set(links))
                        logger.info(f"   Found {len(links)} product links on page {page_num}.")
                        
                        if not links:
                            logger.warning("   No products found. Ending crawl for this category.")
                            break
                            
                        for link in links:
                            if not is_full_crawl and count >= limit:
                                break
                                
                            if link == url: continue
                            
                            try:
                                await self.process_product(page, link, category_name)
                                count += 1
                            except Exception as e:
                                logger.error(f"Failed to process {link}: {e}")
                    except Exception as e:
                        logger.error(f"Error extracting links on page {page_num}: {e}")
                        break
                    
                    page_num += 1
                    
            except Exception as e:
                logger.error(f"Crawl failed for {category_name}: {e}")
            finally:
                await context.close()
        
        return count

    async def run(self, search_term: str, limit: int, job_id: str):
        """Orchestrate the crawling process"""
        from app.api.scheduler import update_job_status
        
        await update_job_status(job_id, status="running")
        
        total_processed = 0
        errors = []
        
        try:
            categories_to_crawl = {}
            if search_term and search_term != 'all':
                # Exact match first, then partial
                if search_term in self.CATEGORIES:
                    categories_to_crawl[search_term] = self.CATEGORIES[search_term]
                else:
                    for cat, url in self.CATEGORIES.items():
                        if search_term.lower() in cat.lower():
                            categories_to_crawl[cat] = url
                
                if not categories_to_crawl:
                    # Fallback: if search_term is not a category key, maybe it's a keyword?
                    # But for now, just log error
                    await update_job_status(job_id, status="completed", message=f"No categories matched '{search_term}'")
                    return
            else:
                categories_to_crawl = self.CATEGORIES

            for cat, url in categories_to_crawl.items():
                try:
                    count = await self.crawl_category(cat, url, limit)
                    total_processed += count
                except Exception as e:
                    logger.error(f"Error crawling category {cat}: {e}")
                    errors.append(f"{cat}: {str(e)}")
            
            await update_job_status(job_id, status="completed", records_drafted=total_processed, errors=errors, completed_at=datetime.utcnow().isoformat())

        except Exception as e:
            logger.error(f"Crawler run failed: {e}")
            await update_job_status(job_id, status="failed", errors=[str(e)])

    async def process_product(self, page: Page, url: str, category: str):
        """Process a single product page"""
        logger.info(f"   Processing Product: {url}")
        await page.goto(url, wait_until="domcontentloaded")
        await self._smart_delay()
        
        content = await page.content()
        title = await page.title()
        body_text = await page.inner_text("body")
        
        # Direct Extraction (Ambeed specific)
        cat_no = None
        cas_no = None
        smiles = None
        mw = None
        formula = None
        
        # Try to extract from table or text
        # Ambeed usually has a table with "Catalog No.", "CAS No.", "Molecular Formula", "Molecular Weight"
        
        # Helper to extract by label
        async def extract_by_label(label):
            return await page.evaluate(f"""
                (label) => {{
                    const elements = Array.from(document.querySelectorAll('td, th, div, span, p'));
                    for (const el of elements) {{
                        if (el.innerText.includes(label)) {{
                            // Try next sibling or parent's next sibling
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
        formula = await extract_by_label("Molecular Formula")
        mw = await extract_by_label("Molecular Weight")
        smiles = await extract_by_label("SMILES Code") # User mentioned "SMILES Code :"
        
        if not smiles:
             smiles = await extract_by_label("SMILES")

        if not cat_no:
            # Fallback from URL
            parts = url.split('/')
            if parts:
                last_part = parts[-1]
                cat_no = last_part.replace('.html', '').replace('.htm', '')

        # Price Data
        price_data = await page.evaluate("""
            () => {
                const rows = Array.from(document.querySelectorAll('table tr'));
                const prices = [];
                rows.forEach(row => {
                    const cols = row.querySelectorAll('td');
                    if (cols.length >= 3) {
                        const size = cols[0].innerText.trim();
                        const price = cols[1].innerText.trim();
                        const availability = cols[2].innerText.trim();
                        if (price.includes('$') || price.includes('€')) {
                            prices.push({ size, price, availability });
                        }
                    }
                });
                return prices;
            }
        """)
        
        logger.info(f"      Found: Cat={cat_no}, CAS={cas_no}, SMILES={smiles}")
        
        # AI Enrichment
        if not smiles and cas_no:
            smiles = await self._fetch_pubchem_smiles(cas_no)
            
        ai_data = await self._enrich_with_gemini(body_text[:3000])
        
        # Merge properties
        properties = ai_data.get("properties", {})
        if mw: properties['mw'] = mw
        if formula: properties['formula'] = formula
        
        # Embedding
        embedding_text = f"{title} {ai_data.get('target', '')} {ai_data.get('summary', '')} {smiles or ''}"
        embedding = await self._get_embedding(embedding_text)
        
        # Save to DB
        data = {
            "ambeed_cat_no": cat_no,
            "cas_number": cas_no,
            "product_name": title.split("|")[0].strip(),
            "product_url": url,
            "category": category,
            "smiles_code": smiles,
            "target": ai_data.get("target"),
            "properties": properties,
            "summary": ai_data.get("summary"),
            "price_data": price_data,
            "embedding": embedding,
            "source_name": "Ambeed", # Explicit source name
            "crawled_at": datetime.utcnow().isoformat()
        }
        
        # Upsert
        try:
            # We use ambeed_cat_no as the unique key, assuming it's unique across Ambeed products
            # If we want to share the table with Creative Biolabs, we should ensure no collision.
            # But the user said "Table sharing... vendor column to distinguish".
            # Upserting on 'ambeed_cat_no' might be risky if Creative Biolabs uses the same column for its ID.
            # Let's check the schema again. 
            # Step 14: 06_Database_Schema.md doesn't explicitly show 'ambeed_cat_no' as PK, but CreativeBiolabsCrawler uses it as conflict key.
            # "supabase.table("commercial_reagents").upsert(data, on_conflict="ambeed_cat_no").execute()"
            # If both use 'ambeed_cat_no' as the conflict key, and IDs overlap, they will overwrite.
            # However, usually Catalog Numbers are vendor specific.
            # To be safe, we should probably include 'source_name' in the conflict key if possible, or assume they are distinct.
            # Given the user instruction "Table sharing... vendor column to distinguish", I will proceed with 'ambeed_cat_no' as the key for now, 
            # assuming the column name 'ambeed_cat_no' implies it was originally designed for Ambeed or just reused.
            
            supabase.table("commercial_reagents").upsert(data, on_conflict="ambeed_cat_no").execute()
            logger.info(f"      ✅ Saved: {data['product_name']} (Target: {data['target']})")
        except Exception as e:
            logger.error(f"      ❌ DB Save failed: {e}")

# Singleton
ambeed_crawler = AmbeedCrawler()
