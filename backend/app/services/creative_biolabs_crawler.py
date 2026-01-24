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
from app.services.ai_refiner import ai_refiner
from app.services.rag_service import rag_service

logger = logging.getLogger(__name__)

class CreativeBiolabsCrawler:
    """Creative Biolabs ADC Reagents Stealth Crawler & AI Enrichment"""
    
    BASE_URL = "https://www.creative-biolabs.com/adc/"
    
    # Categories to crawl (Updated URLs)
    CATEGORIES = {
        "ADC Linker": "https://www.creative-biolabs.com/adc/classify-adc-linker-products-7.htm",
        "ADC Toxin": "https://www.creative-biolabs.com/adc/classify-adc-toxin-products-6.htm",
        "Drug-Linker Complex": "https://www.creative-biolabs.com/adc/classify-drug-linker-complex-products-8.htm"
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
        Analyze the following ADC reagent description and extract structured data.
        
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
        """Generate vector embedding for text using RAG Service (1536 dimensions)"""
        try:
            return await rag_service.generate_embedding(text)
        except Exception as e:
            logger.error(f"Embedding generation failed: {e}")
            return []

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10) -> int:
        """Crawl a specific category with pagination"""
        logger.info(f"🚀 Starting crawl for {category_name}...")
        
        count = 0
        page_num = 1
        max_pages = 116  # Hard limit based on user request
        
        # If limit is very high (e.g. > 1000), we assume full crawl
        is_full_crawl = limit > 1000
        
        async with async_playwright() as p:
            context = await self._init_browser(p)
            page = await context.new_page()
            
            try:
                while True:
                    # Check limits
                    if not is_full_crawl and count >= limit:
                        break
                    if page_num > max_pages:
                        logger.info(f"🏁 Reached max pages ({max_pages}). Stopping.")
                        break

                    # Construct URL with pagination
                    if page_num == 1:
                        url = base_url
                    else:
                        url = f"{base_url}?page={page_num}"
                        
                    logger.info(f"📄 Navigating to Page {page_num}: {url}")
                    
                    try:
                        await page.goto(url, wait_until="networkidle", timeout=60000)
                    except Exception as e:
                        logger.error(f"Page load timeout/error: {e}")
                        break

                    await self._smart_delay()
                    
                    # Stealth: Break every 10 pages
                    if page_num > 1 and page_num % 10 == 0:
                        break_time = 300  # 5 minutes
                        logger.warning(f"☕ Stealth Break! Sleeping for {break_time/60:.1f} minutes...")
                        await asyncio.sleep(break_time)
                    
                    # Log page title for stealth monitoring
                    page_title = await page.title()
                    logger.info(f"   Title: {page_title}")
                    
                    # Extract product links using Regex
                    # Pattern: /adc/[product-name]-[number].htm
                    links = await page.evaluate("""
                        Array.from(document.querySelectorAll('a'))
                            .map(a => a.href)
                            .filter(href => /\\/adc\\/[a-zA-Z0-9-]+-\\d+\\.htm$/.test(href))
                    """)
                    
                    # Remove duplicates and filter out non-product links if any
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
                    
                    # Check for next page (simplified: if we found links, try next page)
                    # In a real scenario, we should check for "Next" button existence.
                    # Given the structure, we can just increment page_num.
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
                    # 'Daily' 모드일 경우 limit를 작게 (예: 10), 'Full Load'일 경우 크게 (예: 1000)
                    # 프론트엔드에서 이미 limit를 보내주므로 그대로 사용
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
        
        # Extract Data
        cat_no = None
        cas_no = None
        
        # Robust Parsing Logic
        lines = body_text.split('\n')
        for line in lines:
            line = line.strip()
            if "Cat. No." in line or "Catalog Number" in line:
                parts = line.split(":")
                if len(parts) > 1:
                    cat_no = parts[-1].strip()
            if "CAS No." in line or "CAS Number" in line:
                parts = line.split(":")
                if len(parts) > 1:
                    cas_no = parts[-1].strip()
                
        if not cat_no:
            # Fallback: try to generate from URL or Title
            # URL format: .../product-name-123.htm -> 123 might be ID, not Cat No.
            # Let's try to find it in title if present (e.g. "Product Name (CAT#: ...)")
            if "(CAT#:" in title:
                try:
                    cat_no = title.split("(CAT#:")[1].split(")")[0].strip()
                except:
                    pass
            
            if not cat_no:
                 cat_no = url.split("/")[-1].replace(".htm", "")
            
        # Extract Price Data (Table parsing)
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
        
        logger.info(f"      Found: Cat={cat_no}, CAS={cas_no}, Prices={len(price_data)}")
        
        # AI Enrichment
        smiles = await self._fetch_pubchem_smiles(cas_no)
        ai_data = await self._enrich_with_gemini(body_text[:3000])
        
        # Embedding
        embedding_text = f"{title} {ai_data.get('target', '')} {ai_data.get('summary', '')} {smiles or ''}"
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
            "summary": ai_data.get("summary"),
            "price_data": price_data,
            "embedding": embedding,
            "source_name": "Creative Biolabs", # Explicit source name
            "crawled_at": datetime.utcnow().isoformat()
        }
        
        # Upsert
        try:
            res = supabase.table("commercial_reagents").upsert(data, on_conflict="ambeed_cat_no").execute()
            logger.info(f"      ✅ Saved: {data['product_name']} (Target: {data['target']})")
            
            # [실시간 연동] 수집 즉시 AI Refiner 호출 (비동기)
            if res.data:
                record_id = res.data[0].get('id')
                if record_id:
                    logger.info(f"      🚀 Triggering Real-time AI Refinement for ID: {record_id}")
                    asyncio.create_task(self._trigger_refinement(res.data[0]))
                    
        except Exception as e:
            logger.error(f"      ❌ DB Save failed: {e}")

    async def _trigger_refinement(self, record: Dict):
        """AI 정제 엔진 비동기 호출"""
        try:
            analysis = await ai_refiner.refine_single_record(record)
            
            if analysis and "error" not in analysis:
                update_data = {
                    "target": analysis.get("target"),
                    "ai_refined": True,
                    "properties": {**record.get("properties", {}), "ai_analysis": analysis}
                }
                supabase.table("commercial_reagents").update(update_data).eq("id", record["id"]).execute()
                logger.info(f"      ✨ Real-time Refinement Success for {record.get('product_name')}")
            else:
                logger.warning(f"      ⚠️ Real-time Refinement failed or skipped for {record.get('product_name')}")
        except Exception as e:
            logger.error(f"      ❌ Real-time Refinement Trigger Error: {e}")

# Singleton
creative_crawler = CreativeBiolabsCrawler()
