import aiohttp
import asyncio
from bs4 import BeautifulSoup
from typing import List, Dict, Any, Optional
from app.core.supabase import supabase
import logging
import re

logger = logging.getLogger(__name__)

class CreativeBiolabsScraper:
    BASE_URL = "https://www.creative-biolabs.com"
    
    async def run(self, search_term: Optional[str] = None, max_pages: int = 3, job_id: Optional[str] = None):
        """
        크롤링 실행 (aiohttp & DB Status)
        """
        from app.api.scheduler import update_job_status
        
        if job_id:
            await update_job_status(job_id, status="running")

        async with aiohttp.ClientSession(headers={"User-Agent": "Mozilla/5.0"}) as session:
            try:
                if search_term:
                    await self.run_search_mode(session, search_term, max_pages, job_id)
                else:
                    await self.run_category_mode(session, max_pages, job_id)
                
                from app.api.scheduler import get_job_from_db
                job = await get_job_from_db(job_id) if job_id else None
                if job and job["status"] != "stopped":
                    await update_job_status(job_id, status="completed", completed_at=datetime.utcnow().isoformat())
            except Exception as e:
                logger.error(f"Crawler Global Error: {e}")
                if job_id:
                    await update_job_status(job_id, status="failed", errors=[str(e)])

    async def run_search_mode(self, session: aiohttp.ClientSession, term: str, max_pages: int, job_id: Optional[str] = None):
        search_url = f"{self.BASE_URL}/adc/search.aspx?q={term}"
        await self.process_list_page(session, search_url, job_id)

    async def run_category_mode(self, session: aiohttp.ClientSession, max_pages: int, job_id: Optional[str] = None):
        categories = [
            "https://www.creative-biolabs.com/adc/products.html",
            "https://www.creative-biolabs.com/adc/adc-cytotoxin.html",
            "https://www.creative-biolabs.com/adc/adc-linkers.html"
        ]
        for cat_url in categories:
            from app.api.scheduler import is_cancelled
            if job_id and await is_cancelled(job_id): break
            await self.process_list_page(session, cat_url, job_id)

    async def process_list_page(self, session: aiohttp.ClientSession, url: str, job_id: Optional[str] = None):
        from app.api.scheduler import update_job_status, get_job_from_db
        try:
            async with session.get(url) as res:
                if res.status != 200: return
                html = await res.text()
                soup = BeautifulSoup(html, "html.parser")
                
                product_links = []
                for a in soup.select("a[href*='/adc/p/'], a[href*='/adc/product/']"):
                    link = a["href"]
                    if not link.startswith("http"):
                        link = self.BASE_URL + link
                    if link not in product_links:
                        product_links.append(link)

                if job_id:
                    job = await get_job_from_db(job_id)
                    current_found = job.get("records_found", 0) if job else 0
                    await update_job_status(job_id, records_found=current_found + len(product_links))

                # 병렬 처리: 5개씩
                for i in range(0, len(product_links), 5):
                    from app.api.scheduler import is_cancelled
                    if job_id and await is_cancelled(job_id): break
                    batch = product_links[i:i+5]
                    tasks = [self.parse_product_page(session, link, job_id) for link in batch]
                    await asyncio.gather(*tasks, return_exceptions=True)
                    
        except Exception as e:
            logger.error(f"List Page Error {url}: {e}")

    async def parse_product_page(self, session: aiohttp.ClientSession, url: str, job_id: Optional[str] = None):
        from app.api.scheduler import update_job_status, get_job_from_db
        try:
            async with session.get(url) as res:
                if res.status != 200: return
                html = await res.text()
                soup = BeautifulSoup(html, "html.parser")
                
                name = soup.select_one("h1").text.strip() if soup.select_one("h1") else "Unknown Product"
                
                # SMILES 보강 (PubChem)
                from app.services.chemical_resolver import chemical_resolver
                smiles = chemical_resolver.fetch_verified_smiles(name)
                
                uniprot_match = re.search(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", html)
                uniprot_id = uniprot_match.group(0) if uniprot_match else None
                
                product_data = {
                    "name": name,
                    "source_url": url,
                    "uniprot_id": uniprot_id,
                    "smiles_code": smiles,
                    "raw_data": {"html_title": soup.title.string if soup.title else None},
                    "status": "available"
                }
                
                supabase.table("commercial_reagents").upsert(product_data, on_conflict="name").execute()
                
                if job_id:
                    job = await get_job_from_db(job_id)
                    current_drafted = job.get("records_drafted", 0) if job else 0
                    await update_job_status(job_id, records_drafted=current_drafted + 1)
                    
        except Exception as e:
            logger.error(f"Product Page Error {url}: {e}")

creative_crawler = CreativeBiolabsScraper()
