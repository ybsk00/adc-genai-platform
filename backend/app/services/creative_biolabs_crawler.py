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

class CreativeBiolabsCrawler:
    """Creative Biolabs ADC Reagents Stealth Crawler & AI Enrichment"""
    
    BASE_URL = "https://www.creative-biolabs.com/adc/"
    
    # Categories to crawl
    CATEGORIES = {
        "ADC Linker-Payload": "https://www.creative-biolabs.com/adc/adc-linker-payloads.htm",
        "ADC Cytotoxin": "https://www.creative-biolabs.com/adc/adc-cytotoxins.htm",
        "ADC Linker": "https://www.creative-biolabs.com/adc/adc-linkers.htm"
    }
    
    def __init__(self):
        self.ua = UserAgent()
        self.request_count = 0
        self.break_threshold = 50  # Break every 50 requests
        
        # Configure Gemini
        genai.configure(api_key=settings.GOOGLE_API_KEY)
        self.model = genai.GenerativeModel('gemini-2.5-flash')

    async def _init_browser(self, p) -> BrowserContext:
        """Initialize browser with stealth settings"""
        user_agent = self.ua.random
        
        browser = await p.chromium.launch(headless=True)
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
            break_time = random.uniform(300, 600)  # 5-10 minutes
            logger.warning(f"☕ Break Time! Sleeping for {break_time/60:.1f} minutes...")
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
        Analyze the following ADC reagent description and extract structured data.
        
        Description: {description[:2000]}
        
        Extract:
        1. Target: The specific antigen target (e.g., HER2, TROP2, CD19). If not specified, return null.
           - Normalize "HER2-directed" to "HER2".
        2. Properties: Extract Purity, Molecular Weight (MW), Formula if present.
        
        Return JSON format:
        {{
            "target": "HER2",
            "properties": {{
                "purity": ">98%",
                "mw": "1234.56",
                "formula": "C50H60N10O12"
            }}
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

    async def crawl_category(self, category_name: str, url: str, limit: int = 10) -> int:
        """Crawl a specific category"""
        logger.info(f"🚀 Starting crawl for {category_name}...")
        
        count = 0
        async with async_playwright() as p:
            context = await self._init_browser(p)
            page = await context.new_page()
            
            try:
                await page.goto(url, wait_until="networkidle")
                await self._smart_delay()
                
                # Get product links (simplified selector, needs adjustment based on actual site)
                # Assuming standard list structure. We might need to adjust this after first run.
                # For now, let's try to find links that look like products.
                # Creative Biolabs usually has a product list table or grid.
                
                # Try to find product links. This is a guess, we might need to refine.
                # Look for links containing 'adc' and ending in '.htm' but not the category page itself
                links = await page.evaluate("""
                    Array.from(document.querySelectorAll('a[href*="/adc/"]'))
                        .map(a => a.href)
                        .filter(href => href.endsWith('.htm') && href.length > 50) 
                """)
                
                # Remove duplicates
                links = list(set(links))
                logger.info(f"Found {len(links)} potential product links.")
                
                for link in links:
                    if count >= limit:
                        break
                        
                    if link == url: continue
                    
                    try:
                        await self.process_product(page, link, category_name)
                        count += 1
                    except Exception as e:
                        logger.error(f"Failed to process {link}: {e}")
                        
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
            # Determine categories to crawl
            categories_to_crawl = {}
            if search_term:
                # Simple case-insensitive partial match for category keys
                for cat, url in self.CATEGORIES.items():
                    if search_term.lower() in cat.lower():
                        categories_to_crawl[cat] = url
                
                if not categories_to_crawl:
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
        logger.info(f"📄 Processing {url}...")
        await page.goto(url, wait_until="domcontentloaded")
        await self._smart_delay()
        
        # Extract Data (Selectors need to be robust)
        # We will extract text content and parse it, as specific selectors might vary.
        
        content = await page.content()
        title = await page.title()
        body_text = await page.inner_text("body")
        
        # Basic extraction using text analysis (since we don't know exact selectors yet)
        # We look for patterns like "Cat. No.:", "CAS No.:"
        
        cat_no = None
        cas_no = None
        
        # Simple parsing logic (can be improved with BeautifulSoup or more specific Playwright selectors)
        lines = body_text.split('\n')
        for line in lines:
            if "Cat. No." in line or "Catalog Number" in line:
                cat_no = line.split(":")[-1].strip()
            if "CAS No." in line or "CAS Number" in line:
                cas_no = line.split(":")[-1].strip()
                
        if not cat_no:
            # Fallback: try to generate from URL or Title
            cat_no = url.split("/")[-1].replace(".htm", "")
            
        logger.info(f"   Found: Cat={cat_no}, CAS={cas_no}")
        
        # AI Enrichment
        smiles = await self._fetch_pubchem_smiles(cas_no)
        ai_data = await self._enrich_with_gemini(body_text[:3000])
        
        # Embedding
        embedding_text = f"{title} {ai_data.get('target', '')} {smiles or ''}"
        embedding = await self._get_embedding(embedding_text)
        
        # Save to DB
        data = {
            "ambeed_cat_no": cat_no,  # Using this column for Catalog No
            "cas_number": cas_no,
            "product_name": title.split("|")[0].strip(),
            "product_url": url,
            "category": category,
            "smiles_code": smiles,
            "target": ai_data.get("target"),
            "properties": ai_data.get("properties"),
            "embedding": embedding,
            "crawled_at": datetime.utcnow().isoformat()
        }
        
        # Upsert
        try:
            supabase.table("commercial_reagents").upsert(data, on_conflict="ambeed_cat_no").execute()
            logger.info(f"   ✅ Saved: {data['product_name']} (Target: {data['target']})")
        except Exception as e:
            logger.error(f"   ❌ DB Save failed: {e}")

# Singleton
creative_crawler = CreativeBiolabsCrawler()
