import httpx
from bs4 import BeautifulSoup
import asyncio
import random
import re
from typing import Optional, List, Dict, Any
from fake_useragent import UserAgent
from app.core.supabase import supabase
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class CreativeBiolabsCrawler:
    BASE_URL = "https://www.creative-biolabs.com"
    # ADC Primary Antibody 카테고리 URL
    CATEGORY_URL = "https://www.creative-biolabs.com/adc/products/primary-antibodies-for-adc-27.htm"
    SEARCH_URL = "https://www.creative-biolabs.com/search"
    
    def __init__(self):
        self.ua = UserAgent()

    def get_headers(self):
        return {
            'User-Agent': self.ua.random,
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
            'Referer': self.BASE_URL,
        }

    async def safe_request(self, client: httpx.AsyncClient, url: str, params: Optional[Dict] = None):
        for i in range(3):
            try:
                await asyncio.sleep(random.uniform(2, 5)) # 매너 딜레이
                headers = self.get_headers()
                res = await client.get(url, headers=headers, params=params)
                if res.status_code == 200:
                    return res
                elif res.status_code in [403, 429]:
                    logger.warning(f"⚠️ Blocked. Cooling down...")
                    await asyncio.sleep(30)
            except Exception as e:
                logger.error(f"Req Error: {e}")
                await asyncio.sleep(5)
        return None

    async def parse_product_page(self, client: httpx.AsyncClient, url: str):
        """상세 페이지 파싱 및 commercial_reagents 저장"""
        res = await self.safe_request(client, url)
        if not res: return

        soup = BeautifulSoup(res.text, 'html.parser')
        
        try:
            h1 = soup.select_one('h1')
            product_name = h1.get_text(strip=True) if h1 else "Unknown Antibody"
            
            properties = {"target": "Unknown"}
            cat_no = None
            supplier = "Creative Biolabs"

            # 상세 스펙 테이블 파싱
            for row in soup.select('table tr'):
                text = row.get_text(" ", strip=True)
                
                if "Catalog No" in text or "Cat No" in text:
                    cat_no = text.split(":")[-1].strip()
                
                if "UniProt" in text:
                    raw_id = text.split("ID")[-1].strip()
                    uniprot_id = re.sub(r'[^A-Z0-9]', '', raw_id) 
                    properties['uniprot_id'] = uniprot_id
                    logger.info(f"🧬 Found UniProt: {uniprot_id}")
                
                if "Target" in text:
                    target = text.split(":")[-1].strip()
                    properties['target'] = target

            if not cat_no:
                cat_no = url.split("/")[-1].replace(".htm", "")

            # 중복 체크 (commercial_reagents 테이블)
            existing = supabase.table("commercial_reagents").select("id").eq("catalog_no", cat_no).execute()
            if existing.data:
                logger.info(f"Skipping duplicate reagent: {cat_no}")
                return

            # commercial_reagents 형식으로 저장
            new_reagent = {
                "name": product_name,
                "category": "antibody",
                "description": f"High affinity antibody targeting {properties['target']} from Creative Biolabs.",
                "supplier": supplier,
                "product_url": url,
                "catalog_no": cat_no,
                "properties": properties
            }
            
            supabase.table("commercial_reagents").insert(new_reagent).execute()
            logger.info(f"✅ Saved to Commercial Reagents: {product_name}")

        except Exception as e:
            logger.error(f"Parse Error {url}: {e}")

    async def run(self, search_term: Optional[str] = None, max_pages: int = 3):
        """
        크롤링 실행
        """
        async with httpx.AsyncClient(timeout=30.0, follow_redirects=True) as client:
            if search_term:
                await self.run_search_mode(client, search_term)
            else:
                await self.run_category_mode(client, max_pages)

    async def run_search_mode(self, client: httpx.AsyncClient, target: str):
        """특정 타겟 검색 모드"""
        logger.info(f"🔎 Search Mode: {target}")
        res = await self.safe_request(client, self.SEARCH_URL, params={'k': f"{target} antibody", 'type': 'product'})
        if res:
            await self.process_list_page(client, res.text)

    async def run_category_mode(self, client: httpx.AsyncClient, max_pages: int):
        """전체 카테고리 순회 모드"""
        logger.info(f"🕸️ Trawl Mode: Scanning ADC antibodies (Max {max_pages} pages)")
        
        for page in range(1, max_pages + 1):
            page_url = f"{self.CATEGORY_URL}?page={page}" if page > 1 else self.CATEGORY_URL
            logger.info(f"📄 Processing Page {page}...")
            
            res = await self.safe_request(client, page_url)
            if not res: break
            
            found_count = await self.process_list_page(client, res.text)
            if found_count == 0:
                logger.info("No more products found. Stopping.")
                break

    async def process_list_page(self, client: httpx.AsyncClient, html_content: str):
        """리스트 페이지에서 제품 링크 추출 후 방문"""
        soup = BeautifulSoup(html_content, 'html.parser')
        links = []
        for a in soup.select('.pro_list_title a'):
            href = a.get('href')
            if href:
                full_url = self.BASE_URL + href if href.startswith('/') else href
                links.append(full_url)
        
        logger.info(f"Found {len(links)} products on this page.")
        for link in links:
            await self.parse_product_page(client, link)
            
        return len(links)

creative_crawler = CreativeBiolabsCrawler()
