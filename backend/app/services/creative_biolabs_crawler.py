"""
Creative Biolabs ADC 제품 크롤러
카테고리 페이지 및 키워드 검색 기반 크롤링
"""
import aiohttp
import asyncio
from datetime import datetime
from bs4 import BeautifulSoup
from typing import List, Dict, Any, Optional
from app.core.supabase import supabase
import logging
import re

logger = logging.getLogger(__name__)

class CreativeBiolabsScraper:
    BASE_URL = "https://www.creative-biolabs.com"
    
    # 실제 제품 카테고리 URL (사이트 구조 분석 결과)
    CATEGORY_URLS = [
        "https://www.creative-biolabs.com/adc/target-list-customized-adcs-1.htm",
        "https://www.creative-biolabs.com/adc/classify-druglnk-products-2.htm",
        "https://www.creative-biolabs.com/adc/target-list-anti-ab-adcs-3.htm",
        "https://www.creative-biolabs.com/adc/classify-anti-drug-abs-4.htm",
        "https://www.creative-biolabs.com/adc/target-list-customized-fluoroab-9.htm",
        "https://www.creative-biolabs.com/adc/classify-fluorescent-dyes-10.htm"
    ]
    
    # ADC 관련 검색 키워드
    SEARCH_TERMS = [
        "HER2", "TROP2", "EGFR", "CD19", "CD22", "CD33", 
        "MMAE", "DM1", "DXd", "vedotin", "trastuzumab"
    ]
    
    async def run(self, search_term: Optional[str] = None, max_pages: int = 5, job_id: Optional[str] = None):
        """
        크롤링 실행 (aiohttp & DB Status)
        """
        from app.api.scheduler import update_job_status, get_job_from_db
        
        if job_id:
            await update_job_status(job_id, status="running")

        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.5"
        }

        async with aiohttp.ClientSession(headers=headers) as session:
            try:
                if search_term:
                    # 특정 키워드로 검색
                    await self.run_search_mode(session, search_term, job_id)
                else:
                    # 카테고리 페이지 + 키워드 검색 모드
                    await self.run_full_crawl_mode(session, job_id)
                
                job = await get_job_from_db(job_id) if job_id else None
                if job and job["status"] != "stopped":
                    await update_job_status(job_id, status="completed", completed_at=datetime.utcnow().isoformat())
            except Exception as e:
                logger.error(f"Crawler Global Error: {e}")
                if job_id:
                    await update_job_status(job_id, status="failed", errors=[str(e)])

    async def run_search_mode(self, session: aiohttp.ClientSession, term: str, job_id: Optional[str] = None):
        """키워드 검색 모드"""
        search_url = f"{self.BASE_URL}/adc/search.aspx?q={term}"
        logger.info(f"🔍 Searching for: {term}")
        await self.process_list_page(session, search_url, job_id)

    async def run_full_crawl_mode(self, session: aiohttp.ClientSession, job_id: Optional[str] = None):
        """전체 크롤링 모드: 카테고리 페이지 + 키워드 검색"""
        from app.api.scheduler import is_cancelled, update_job_status
        
        # 1. 카테고리 페이지 크롤링
        logger.info("📂 Crawling category pages...")
        for cat_url in self.CATEGORY_URLS:
            if job_id and await is_cancelled(job_id):
                break
            await self.process_list_page(session, cat_url, job_id)
            await asyncio.sleep(1)  # Rate limiting
        
        # 2. 키워드 검색 크롤링
        logger.info("🔍 Crawling by search keywords...")
        for term in self.SEARCH_TERMS:
            if job_id and await is_cancelled(job_id):
                break
            search_url = f"{self.BASE_URL}/adc/search.aspx?q={term}"
            await self.process_list_page(session, search_url, job_id)
            await asyncio.sleep(1)  # Rate limiting

    async def process_list_page(self, session: aiohttp.ClientSession, url: str, job_id: Optional[str] = None):
        """목록 페이지에서 제품 링크 추출 및 처리"""
        from app.api.scheduler import update_job_status, get_job_from_db, is_cancelled
        
        try:
            async with session.get(url, timeout=aiohttp.ClientTimeout(total=30)) as res:
                if res.status != 200:
                    logger.warning(f"Page load failed: {url} -> {res.status}")
                    return
                    
                html = await res.text()
                soup = BeautifulSoup(html, "html.parser")
                
                # 제품 링크 추출 (다양한 패턴)
                product_links = set()
                
                # 패턴 1: 제품 상세 페이지 링크
                for a in soup.select("a[href*='.htm']"):
                    href = a.get('href', '')
                    if '/adc/' in href and not any(x in href for x in ['search', 'home', 'contact', 'about']):
                        if not href.startswith("http"):
                            href = self.BASE_URL + href
                        product_links.add(href)
                
                # 패턴 2: 제품 카드/리스트 아이템
                for a in soup.select(".product-item a, .item-box a, .product-list a"):
                    href = a.get('href', '')
                    if href and not href.startswith("http"):
                        href = self.BASE_URL + href
                    if href:
                        product_links.add(href)

                logger.info(f"Found {len(product_links)} product links on {url[:50]}...")

                if job_id:
                    job = await get_job_from_db(job_id)
                    current_found = job.get("records_found", 0) if job else 0
                    await update_job_status(job_id, records_found=current_found + len(product_links))

                # 병렬 처리: 5개씩
                product_list = list(product_links)
                for i in range(0, len(product_list), 5):
                    if job_id and await is_cancelled(job_id):
                        break
                    batch = product_list[i:i+5]
                    tasks = [self.parse_product_page(session, link, job_id) for link in batch]
                    await asyncio.gather(*tasks, return_exceptions=True)
                    await asyncio.sleep(0.5)  # Rate limiting
                    
        except Exception as e:
            logger.error(f"List Page Error {url}: {e}")

    async def parse_product_page(self, session: aiohttp.ClientSession, url: str, job_id: Optional[str] = None):
        """제품 상세 페이지 파싱 및 저장"""
        from app.api.scheduler import update_job_status, get_job_from_db
        
        try:
            async with session.get(url, timeout=aiohttp.ClientTimeout(total=20)) as res:
                if res.status != 200:
                    return
                    
                html = await res.text()
                soup = BeautifulSoup(html, "html.parser")
                
                # 제품명 추출
                name = "Unknown Product"
                for selector in ["h1", ".product-name", ".product-title", "title"]:
                    elem = soup.select_one(selector)
                    if elem:
                        name = elem.text.strip()
                        break
                
                # 제품 설명 추출
                description = ""
                for selector in [".product-description", ".description", ".content", "p"]:
                    elem = soup.select_one(selector)
                    if elem:
                        description = elem.text.strip()[:500]
                        break
                
                # 카탈로그 번호 추출
                catalog_match = re.search(r"(CAT|CB|BCA|ADC)[A-Z0-9\-]+", html)
                catalog_number = catalog_match.group(0) if catalog_match else None
                
                # SMILES 보강 (PubChem)
                smiles = None
                try:
                    from app.services.chemical_resolver import chemical_resolver
                    smiles = chemical_resolver.fetch_verified_smiles(name)
                except Exception:
                    pass
                
                # UniProt ID 추출
                uniprot_match = re.search(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}", html)
                uniprot_id = uniprot_match.group(0) if uniprot_match else None
                
                product_data = {
                    "name": name[:200],
                    "source_url": url,
                    "catalog_number": catalog_number,
                    "description": description,
                    "uniprot_id": uniprot_id,
                    "smiles_code": smiles,
                    "raw_data": {"html_title": soup.title.string if soup.title else None},
                    "source": "creative_biolabs",
                    "status": "available"
                }
                
                # commercial_reagents 테이블에 저장 (upsert)
                supabase.table("commercial_reagents").upsert(product_data, on_conflict="name").execute()
                
                if job_id:
                    job = await get_job_from_db(job_id)
                    current_drafted = job.get("records_drafted", 0) if job else 0
                    await update_job_status(job_id, records_drafted=current_drafted + 1)
                
                logger.info(f"✅ Saved: {name[:50]}...")
                    
        except Exception as e:
            logger.error(f"Product Page Error {url}: {e}")

# 싱글톤 인스턴스
creative_crawler = CreativeBiolabsScraper()
