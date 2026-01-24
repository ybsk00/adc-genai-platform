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
        
        # [고도화] InChI 필드에서도 SMILES 추출 시도 (비용 절감)
        if not smiles:
            inchi = await extract_by_label("InChI")
            if inchi:
                logger.info(f"      Found InChI: {inchi[:30]}... (Will attempt to use this)")
                # InChI가 있으면 나중에 Chemical Resolver가 처리할 수 있도록 저장
        
        # [고도화] 페이지 하단 텍스트 영역에서 SMILES 패턴 직접 검색
        if not smiles:
            smiles_from_text = await page.evaluate("""
                () => {
                    const bodyText = document.body.innerText;
                    // SMILES 패턴 매칭 (단순화된 정규식)
                    const match = bodyText.match(/SMILES Code\\s*:\\s*([A-Za-z0-9@#\\(\\)\\[\\]\\/\\\\=+-]+)/);
                    return match ? match[1].trim() : null;
                }
            """)
            if smiles_from_text:
                smiles = smiles_from_text
                logger.info(f"      ✅ Found SMILES via Text Pattern: {smiles[:30]}...")

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
            res = supabase.table("commercial_reagents").upsert(data, on_conflict="ambeed_cat_no").execute()
            logger.info(f"      ✅ Saved: {data['product_name']} (Target: {data['target']})")
            
            # [실시간 연동] 수집 즉시 AI Refiner 호출 (비동기)
            if res.data:
                record_id = res.data[0].get('id')
                if record_id:
                    logger.info(f"      🚀 Triggering Real-time AI Refinement for ID: {record_id}")
                    # 정제 작업 때문에 수집 속도가 느려지지 않게 비동기로 처리
                    asyncio.create_task(self._trigger_refinement(res.data[0]))
                    
        except Exception as e:
            logger.error(f"      ❌ DB Save failed: {e}")

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
