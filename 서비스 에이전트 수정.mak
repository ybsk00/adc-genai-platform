crawler_stealth.py (수정 필수)
이유: 현재 코드는 "가격표"를 가져오는 로직이 Placeholder로 비어 있습니다. 아까 보내주신 스크린샷(image_ec94e5.png)의 구조(SMILES Code, Price)를 정확히 긁어오도록 파싱 로직을 정교화해야 합니다.

덮어쓰기 코드:

Python
import requests
from bs4 import BeautifulSoup
import time
import random
from fake_useragent import UserAgent
from supabase import create_client, Client
import os
from dotenv import load_dotenv
import logging
import json

# Load environment variables
load_dotenv()

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class AmbeedStealthScraper:
    BASE_URL = "https://www.ambeed.com"
    
    def __init__(self):
        self.ua = UserAgent()
        self.session = requests.Session()
        
        # Initialize Supabase client
        url: str = os.environ.get("SUPABASE_URL")
        key: str = os.environ.get("SUPABASE_SERVICE_KEY")
        if not url or not key:
            raise ValueError("Supabase credentials not found in environment variables")
        self.supabase: Client = create_client(url, key)
        
    def get_headers(self, referer=None):
        headers = {
            'User-Agent': self.ua.random,
            'Accept-Language': 'en-US,en;q=0.9',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
            'Connection': 'keep-alive',
            'Upgrade-Insecure-Requests': '1',
        }
        if referer:
            headers['Referer'] = referer
        return headers

    def safe_request(self, url, referer=None):
        """
        Anti-blocking request wrapper
        """
        retries = 3
        backoff = 30
        
        for i in range(retries):
            try:
                # Human-like random delay
                sleep_time = random.uniform(2, 5)
                time.sleep(sleep_time)
                
                logger.info(f"🕵️ Scraping: {url} (Delay: {sleep_time:.2f}s)")
                
                response = self.session.get(
                    url, 
                    headers=self.get_headers(referer), 
                    timeout=20
                )
                
                if response.status_code in [403, 429, 503]:
                    logger.warning(f"⚠️ Blocked ({response.status_code}). Cooling down for {backoff}s...")
                    time.sleep(backoff)
                    backoff *= 2
                    self.session = requests.Session()
                    continue
                
                response.raise_for_status()
                return response
            
            except Exception as e:
                logger.error(f"❌ Error fetching {url}: {str(e)}")
                time.sleep(5)
                
        return None

    def parse_detail_page(self, url, category, referer_url):
        """
        [Updated] Parse product detail page strictly based on Ambeed layout
        """
        res = self.safe_request(url, referer=referer_url)
        if not res: return
        
        soup = BeautifulSoup(res.text, 'html.parser')
        
        details = {
            'category': category,
            'product_url': url,
            'price_data': {} 
        }
        
        try:
            # 1. Product Name
            title_elem = soup.select_one('h1')
            if title_elem:
                details['product_name'] = title_elem.get_text(strip=True)

            # 2. Specifications Table (스크린샷 기반 파싱)
            # Ambeed는 보통 <th>Label</th><td>Value</td> 또는 div 구조를 씀.
            # 텍스트 기반으로 찾는 것이 가장 안전함.
            
            # 전체 텍스트에서 키워드 추출 시도
            for row in soup.find_all(['tr', 'li', 'div']): 
                text = row.get_text(" ", strip=True)
                
                if "Catalog No" in text and ":" in text:
                    details['ambeed_cat_no'] = text.split(":")[-1].strip()
                elif "CAS No" in text and ":" in text:
                    details['cas_number'] = text.split(":")[-1].strip()
                elif "Formula" in text and ":" in text:
                    details['formula'] = text.split(":")[-1].strip()
                elif "M.W" in text and ":" in text:
                    details['molecular_weight'] = text.split(":")[-1].strip()
                elif "SMILES" in text and ":" in text:
                    # SMILES 코드는 매우 중요
                    raw_smiles = text.split(":")[-1].strip()
                    details['smiles_code'] = raw_smiles

            # 3. Price Parsing (가격표)
            # 보통 table class="price-table" 같은 형식이므로, '$'가 포함된 행을 찾음
            price_list = []
            for row in soup.find_all('tr'):
                row_text = row.get_text(strip=True)
                if "$" in row_text and ("mg" in row_text or "g" in row_text):
                    # 예: "1mg $100" -> 간단히 텍스트 전체 저장 (나중에 AI가 해석)
                    price_list.append(row_text)
            
            if price_list:
                details['price_data'] = {"raw_pricing": price_list}

            # 4. Save if valid
            if 'ambeed_cat_no' in details:
                self.save_to_db(details)
            else:
                logger.warning(f"Skipping {url}: Could not find Catalog No")

        except Exception as e:
            logger.error(f"Error parsing {url}: {str(e)}")

    def save_to_db(self, data):
        """Upsert to Supabase commercial_reagents table"""
        try:
            # Check existing
            existing = self.supabase.table('commercial_reagents')\
                .select('id')\
                .eq('ambeed_cat_no', data['ambeed_cat_no'])\
                .execute()
                
            if existing.data:
                self.supabase.table('commercial_reagents')\
                    .update(data)\
                    .eq('ambeed_cat_no', data['ambeed_cat_no'])\
                    .execute()
                logger.info(f"🔄 Updated {data.get('product_name')} ({data['ambeed_cat_no']})")
            else:
                self.supabase.table('commercial_reagents')\
                    .insert(data)\
                    .execute()
                logger.info(f"✅ Inserted {data.get('product_name')} ({data['ambeed_cat_no']})")

        except Exception as e:
            logger.error(f"DB Error: {str(e)}")

    def run(self, category="Payload"):
        # ... (기존 run 로직 유지) ...
        pass
2. 🔍 rag_service.py (수정 필수)
이유: Commercial Agent가 "혹시 MMAE 접합체 있어요?"라고 물어봤을 때, DB를 검색해주는 기능(search_commercial)이 빠져 있습니다. 이걸 추가해야 에이전트가 데이터를 찾을 수 있습니다.

덮어쓰기 코드:

Python
import os
from typing import List, Dict, Any
from openai import AsyncOpenAI
from app.core.supabase import supabase

class RAGService:
    def __init__(self):
        self.api_key = os.getenv("OPENAI_API_KEY")
        self.client = AsyncOpenAI(api_key=self.api_key)

    async def generate_embedding(self, text: str) -> List[float]:
        """Generate embedding for text using OpenAI"""
        response = await self.client.embeddings.create(
            input=text,
            model="text-embedding-3-small"
        )
        return response.data[0].embedding

    # [NEW] 상용 시약 검색 기능 (Commercial Agent가 사용)
    async def search_commercial(self, query: str, limit: int = 3) -> List[Dict[str, Any]]:
        """
        Search for commercial reagents (Ambeed Data).
        First try: Text Match (Product Name / Cat No)
        Second try: Semantic Search (Description/Category)
        """
        try:
            # 1. 텍스트 검색 (정확도 우선)
            # Supabase의 textSearch 기능 활용 (product_name 컬럼)
            text_result = supabase.table("commercial_reagents")\
                .select("*")\
                .textSearch("product_name", query, config="english")\
                .limit(limit)\
                .execute()
            
            if text_result.data:
                return text_result.data

            # 2. 결과가 없으면 의미 기반 검색 (Semantic Search)
            # (만약 commercial_reagents 테이블에 embedding 컬럼을 만들었다면 아래 로직 사용)
            # embedding = await self.generate_embedding(query)
            # rpc_result = supabase.rpc("match_commercial", { ... }).execute()
            # return rpc_result.data
            
            return []
            
        except Exception as e:
            print(f"Error searching commercial DB: {str(e)}")
            return []

    # ... (기존 index_golden_set_item, search 메서드 유지) ...
    async def index_golden_set_item(self, item_id: str, description: str, properties: Dict[str, Any]):
        # (기존 코드 유지)
        pass

    async def search(self, query: str, top_k: int = 5) -> List[Dict[str, Any]]:
        # (기존 코드 유지)
        pass

rag_service = RAGService()