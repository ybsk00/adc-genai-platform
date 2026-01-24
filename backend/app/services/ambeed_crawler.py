import asyncio
import random
import logging
import json
import time
from typing import Dict, List, Optional, Any
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

class AmbeedCrawler:
    """Ambeed ADC Reagents Stealth Crawler & AI Enrichment"""
    
    BASE_URL = "https://www.ambeed.com"
    
    # Categories to crawl (Direct URLs)
    CATEGORIES = {
        "ADC Toxins": "https://www.ambeed.com/adc-toxins.html",
        "Conjugate": "https://www.ambeed.com/antibody-drug-conjugate-adc-related.html",
        "Payload-Linker": "https://www.ambeed.com/payload-linker-conjugate.html",
    }
    
    def __init__(self):
        self.ua = UserAgent()
        self.request_count = 0
        self.break_threshold = 100  # Break every 100 requests
        self.batch_size = 5  # Parallel processing batch size
        
        # Configure Gemini
        # Updated to gemini-2.5-flash as per instruction
        genai.configure(api_key=settings.GOOGLE_API_KEY)
        self.model_id = 'gemini-2.5-flash'
        try:
             self.model = genai.GenerativeModel(self.model_id)
        except Exception:
             logger.warning(f"⚠️ Model {self.model_id} not found, falling back to gemini-1.5-flash")
             self.model = genai.GenerativeModel('gemini-1.5-flash')

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

    async def _enrich_with_gemini(self, description: str, current_data: Dict) -> Dict:
        """Extract Target and Attributes using Gemini 2.5 Flash"""
        
        # 4. Conditional Skip Logic: Cost Savings
        # If we have all core physical properties and description is short, skip AI
        has_core_info = all([
            current_data.get('cas_number'),
            current_data.get('properties', {}).get('mw'),
            current_data.get('properties', {}).get('formula'),
            current_data.get('smiles_code')
        ])
        
        # If text is too short to contain hidden info and we have core info, skip
        if has_core_info and len(description) < 200:
            logger.info("⏩ Skipping AI Enrichment (Core info present & text short)")
            return {"target": None, "properties": {}, "summary": "Auto-generated from structured data."}

        # 3. Enhanced Target Normalization Prompt
        prompt = f"""
        Analyze the following Ambeed reagent description and extract structured data.
        
        Description: {description[:3000]}
        
        Task:
        1. **Target Normalization**: Identify the specific antigen target (e.g., HER2, TROP2, CD19). 
           - **CRITICAL**: Normalize synonyms to standard gene symbols.
           - Examples: "ERBB2" -> "HER2", "CD340" -> "HER2", "TNFSF11" -> "RANKL".
           - If not specified, return null.
        2. **Properties**: Extract Purity, Molecular Weight (MW), Formula if present in text but missed in parsing.
        3. **Summary**: A concise 3-sentence summary suitable for ADC researchers.
        
        Output JSON:
        {{
            "target": "HER2",
            "properties": {{
                "purity": ">98%",
                "mw": "1234.56",
                "formula": "C50H..."
            }},
            "summary": "..."
        }}
        """
        
        try:
            # 1. Update to JSON mode for robustness
            generation_config = genai.GenerationConfig(
                response_mime_type="application/json"
            )
            
            response = await self.model.generate_content_async(
                prompt,
                generation_config=generation_config
            )
            
            text = response.text.strip()
            return json.loads(text)
        except Exception as e:
            logger.error(f"Gemini enrichment failed: {e}")
            return {"target": None, "properties": {}}

    async def _get_embedding(self, text: str) -> Optional[List[float]]:
        """Generate vector embedding for text using RAG Service (768 dimensions)"""
        try:
            embedding = await rag_service.generate_embedding(text)
            if not embedding:
                return None
            if len(embedding) != 768:
                logger.warning(f"⚠️ Generated embedding has {len(embedding)} dimensions, expected 768.")
            return embedding
        except Exception as e:
            logger.error(f"Embedding generation failed: {e}")
            return None

    async def crawl_category(self, category_name: str, base_url: str, limit: int = 10) -> int:
        """Crawl a specific category with pagination and batch processing"""
        logger.info(f"🚀 Starting crawl for {category_name}...")
        
        count = 0
        page_num = 1
        max_pages = 116  # Hard limit
        
        is_full_crawl = limit > 1000
        batch_buffer = [] # Buffer for batch processing
        
        async with async_playwright() as p:
            context = await self._init_browser(p)
            page = await context.new_page()
            
            try:
                while True:
                    if not is_full_crawl and count >= limit:
                        break
                    if page_num > max_pages:
                        break

                    # Ambeed pagination
                    url = base_url if page_num == 1 else f"{base_url}&page={page_num}"
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
                                # 1. Extract HTML Data Only (Fast)
                                raw_data = await self._extract_html_data(page, link, category_name)
                                if raw_data:
                                    batch_buffer.append(raw_data)
                                    count += 1
                                
                                # 2. Process Batch if full
                                if len(batch_buffer) >= self.batch_size:
                                    logger.info(f"⚡ Processing Batch of {len(batch_buffer)} items...")
                                    await self._process_batch_enrichment(batch_buffer)
                                    batch_buffer = [] # Clear buffer
                                    
                            except Exception as e:
                                logger.error(f"Failed to extract {link}: {e}")
                        
                        # Process remaining items in buffer at end of page
                        if batch_buffer:
                            await self._process_batch_enrichment(batch_buffer)
                            batch_buffer = []

                    except Exception as e:
                        logger.error(f"Error extracting links on page {page_num}: {e}")
                        break
                    
                    page_num += 1
                    
            except Exception as e:
                logger.error(f"Crawl failed for {category_name}: {e}")
            finally:
                # Flush any remaining items
                if batch_buffer:
                    await self._process_batch_enrichment(batch_buffer)
                await context.close()
        
        return count

    async def _extract_html_data(self, page: Page, url: str, category: str) -> Optional[Dict]:
        """Navigate to product page and extract raw HTML data"""
        logger.info(f"   extracting HTML: {url}")
        await page.goto(url, wait_until="domcontentloaded")
        await self._smart_delay()
        
        title = await page.title()
        body_text = await page.inner_text("body")
        
        # Helper to extract by label
        async def extract_by_label(label):
            return await page.evaluate(f"""
                (label) => {{
                    const elements = Array.from(document.querySelectorAll('td, th, div, span, p, b, strong, li'));
                    for (const el of elements) {{
                        if (el.children.length > 2 || el.innerText.length > 300) continue;
                        const text = el.innerText.trim();
                        if (text.includes(label)) {{
                            const safeLabel = label.replace(/[.*+?^${{}}()|[\\]\\\\]/g, '\\\\$&');
                            const regex = new RegExp(safeLabel + '\\\\s*[:\\\\uff1a]\\\\s*([^\\\\n\\\\r]+)', 'i');
                            const match = text.match(regex);
                            if (match) return match[1].trim();
                            if (text.replace(':', '').trim() === label || text.length < label.length + 5) {{
                                if (el.nextElementSibling) return el.nextElementSibling.innerText.trim();
                                if (el.parentElement && el.parentElement.nextElementSibling) {{
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
        if not smiles: smiles = await extract_by_label("SMILES")
        
        if not smiles:
            smiles_from_text = await page.evaluate("""
                () => {
                    const bodyText = document.body.innerText;
                    const match = bodyText.match(/SMILES Code\\s*:\\s*([^\\n]+)/);
                    if (match) return match[1].trim();
                    return null;
                }
            """)
            if smiles_from_text: smiles = smiles_from_text

        if not cat_no:
            parts = url.split('/')
            if parts: cat_no = parts[-1].replace('.html', '').replace('.htm', '')

        price_data = await page.evaluate("""
            () => {
                const rows = Array.from(document.querySelectorAll('table tr'));
                const prices = [];
                rows.forEach(row => {
                    const cols = row.querySelectorAll('td');
                    if (cols.length >= 3) {
                        const size = cols[0].innerText.trim();
                        const price = cols[1].innerText.trim();
                        if (price.includes('$') || price.includes('€')) {
                            prices.push({ size, price, availability: cols[2].innerText.trim() });
                        }
                    }
                });
                return prices;
            }
        """)
        
        return {
            "ambeed_cat_no": cat_no,
            "cas_number": cas_no,
            "product_name": title.split("|")[0].strip(),
            "product_url": url,
            "category": category,
            "smiles_code": smiles,
            "properties": {
                "formula": formula,
                "mw": mw
            },
            "price_data": price_data,
            "body_text": body_text, # Passed for AI analysis
            "source_name": "Ambeed",
            "crawled_at": datetime.utcnow().isoformat()
        }

    async def _process_batch_enrichment(self, batch_data: List[Dict]):
        """Process a batch of items with AI in parallel"""
        tasks = [self._enrich_and_save_single(item) for item in batch_data]
        await asyncio.gather(*tasks)

    async def _enrich_and_save_single(self, raw_data: Dict):
        """Enrich a single item and save to DB"""
        try:
            body_text = raw_data.pop("body_text", "")
            
            # PubChem Fallback
            if not raw_data['smiles_code'] and raw_data['cas_number']:
                raw_data['smiles_code'] = await self._fetch_pubchem_smiles(raw_data['cas_number'])

            # AI Enrichment (with conditional skip)
            ai_data = await self._enrich_with_gemini(body_text, raw_data)
            
            # Merge AI data
            raw_data["target"] = ai_data.get("target")
            raw_data["summary"] = ai_data.get("summary")
            
            ai_props = ai_data.get("properties", {})
            if ai_props.get("mw"): raw_data["properties"]["mw"] = ai_props["mw"]
            if ai_props.get("formula"): raw_data["properties"]["formula"] = ai_props["formula"]
            if ai_props.get("purity"): raw_data["properties"]["purity"] = ai_props["purity"]

            # Generate Embedding
            embedding_text = f"{raw_data['product_name']} {raw_data.get('target', '')} {raw_data.get('summary', '')} {raw_data.get('smiles_code', '')}"
            raw_data["embedding"] = await self._get_embedding(embedding_text)

            # Save to DB
            res = supabase.table("commercial_reagents").upsert(raw_data, on_conflict="ambeed_cat_no").execute()
            logger.info(f"      ✅ Saved: {raw_data['product_name']}")
            
            # Trigger Refinement
            if res.data:
                logger.info(f"      🚀 Triggering Real-time AI Refinement for ID: {res.data[0].get('id')}")
                asyncio.create_task(self._trigger_refinement(res.data[0]))
                
        except Exception as e:
            logger.error(f"      ❌ Processing failed for {raw_data.get('product_name')}: {e}")

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

    async def _trigger_refinement(self, record: Dict):
        """AI 정제 엔진 비동기 호출"""
        try:
            # commercial_reagents 테이블의 데이터를 golden_set_library 형식으로 변환하여 전달하거나
            # ai_refiner가 commercial_reagents도 지원하도록 확장 필요.
            # 현재 ai_refiner.py는 golden_set_library를 기준으로 작성되어 있음.
            # 사장님 지시는 "데이터가 DB에 인서트되는 즉시 ai_refined 상태가 False에서 True로 변하며 분류가 완료되어야 함"
            # commercial_reagents 테이블에도 ai_refined 컬럼이 있는지 확인 필요.
            # 만약 없다면 golden_set_library로 데이터를 넘겨주는 로직이 필요할 수 있음.
            
            # 일단 ai_refiner의 refine_single_record를 호출 시도
            # (ai_refiner.py 83행: refine_single_record(record))
            analysis = await ai_refiner.refine_single_record(record)
            
            if analysis and "error" not in analysis:
                # 결과 업데이트 (commercial_reagents 테이블 기준)
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
ambeed_crawler = AmbeedCrawler()
