# 🕵️‍♂️ Project Treasure Trove: Ambeed Stealth Crawler Spec

**Objective**: Ambeed.com의 ADC 관련 제품(Payload, Linker, Conjugate) 전량을 크롤링합니다.
**Critical Requirement**: 웹사이트의 보안 차단(IP Block/403 Error)을 우회하기 위해 User-Agent 로테이션, 랜덤 딜레이, 세션 유지 등의 스텔스 기술을 반드시 적용해야 합니다.

**Target Hub**: Ambeed ADC Related Products

## 1. 🛡️ Anti-Blocking Strategy (우회 기술 지침)
**개발자 필독**: 단순 `requests.get()` 루프를 돌리면 100% 차단당합니다. 아래 기술을 모두 적용하십시오.

### Random User-Agent Rotation
매 요청마다 브라우저 정보(User-Agent)를 최신 Chrome, Firefox, Safari 등으로 랜덤하게 변경하여 "여러 사람이 접속하는 척" 해야 합니다.
*   **라이브러리 추천**: `fake-useragent`

### Human-like Delays (Random Sleep)
페이지 이동 시 기계적으로 1초마다 요청하지 말고, **2초~5초 사이의 랜덤한 시간**을 쉬어야 합니다.

### Referer Header Manipulation
상세 페이지로 들어갈 때, 반드시 **"이전 목록 페이지에서 클릭해서 들어왔음"**을 증명하는 `Referer` 헤더를 심어야 합니다.

### Error Handling & Backoff
403(Forbidden)이나 429(Too Many Requests) 에러 발생 시, 스크립트를 멈추지 말고 **30초 이상 대기 후 재시도(Exponential Backoff)** 하십시오.

## 2. 💾 Database Schema (PostgreSQL)
```sql
CREATE TABLE IF NOT EXISTS commercial_reagents (
    id uuid DEFAULT gen_random_uuid() PRIMARY KEY,
    ambeed_cat_no text UNIQUE NOT NULL,
    product_name text,
    category text, -- 'Payload', 'Linker', 'Conjugate'
    cas_number text,
    smiles_code text, -- 상세 페이지 'Product Details'에서 파싱
    molecular_weight text,
    formula text,
    price_data jsonb,
    stock_status text,
    product_url text,
    crawled_at timestamptz DEFAULT now(),
    embedding vector(1536) 
);
```

## 3. 🐍 Stealth Crawler Python Logic (services/crawler_stealth.py)
아래 코드는 차단을 회피하기 위한 핵심 로직이 포함된 템플릿입니다.

```python
import requests
from bs4 import BeautifulSoup
import time
import random
from fake_useragent import UserAgent # pip install fake-useragent

class AmbeedStealthScraper:
    BASE_URL = "https://www.ambeed.com"
    
    def __init__(self):
        self.ua = UserAgent()
        self.session = requests.Session() # 세션 유지 (쿠키 관리)
        
    def get_headers(self, referer=None):
        headers = {
            'User-Agent': self.ua.random, # 매번 다른 브라우저인 척
            'Accept-Language': 'en-US,en;q=0.9',
            'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        }
        if referer:
            headers['Referer'] = referer
        return headers

    def safe_request(self, url, referer=None):
        """차단 방지를 위한 안전한 요청 함수"""
        retries = 3
        for i in range(retries):
            try:
                # 1. 사람처럼 랜덤하게 쉬기 (2~5초)
                sleep_time = random.uniform(2, 5)
                time.sleep(sleep_time)
                
                print(f"🕵️ Scraping: {url} (Delay: {sleep_time:.2f}s)")
                
                response = self.session.get(
                    url, 
                    headers=self.get_headers(referer), 
                    timeout=15
                )
                
                # 2. 차단 감지 시 대기
                if response.status_code in [403, 429, 503]:
                    print(f"⚠️ Blocked ({response.status_code}). Cooling down for 60s...")
                    time.sleep(60) # 1분 휴식
                    self.session = requests.Session() # 세션 초기화 (새 신분 세탁)
                    continue
                
                response.raise_for_status()
                return response
            
            except Exception as e:
                print(f"❌ Error: {e}")
                time.sleep(5)
        return None

    def parse_detail_page(self, url, category, referer_url):
        res = self.safe_request(url, referer=referer_url)
        if not res: return
        
        soup = BeautifulSoup(res.text, 'html.parser')
        
        # [핵심] 상세 스펙 테이블 파싱 (스크린샷 기반)
        details = {}
        # 실제 사이트 구조에 맞춰 Selector 수정 필요
        # 예: <td class="label">SMILES Code :</td><td class="value">...</td>
        for row in soup.select('tr'): 
            text = row.get_text(strip=True)
            if "SMILES Code" in text:
                details['smiles'] = text.split(':')[-1].strip()
            elif "CAS No." in text:
                details['cas'] = text.split(':')[-1].strip()
        
        # DB 저장 로직 호출
        self.save_to_db(details)

    def run(self):
        # ... (카테고리 순회 로직은 이전과 동일, safe_request 사용 필수) ...
        pass
```

## 4. 🧠 RAG Integration (활용 시나리오)
크롤링 된 데이터는 단순 저장이 아니라 벡터 임베딩 되어야 합니다.

*   **Trigger**: 크롤링이 완료된 직후 (`save_to_db` 이후).
*   **Action**: SMILES Code와 Product Description을 텍스트로 합쳐서 OpenAI Embedding API 호출.
*   **Query**: 사용자가 "MMAE랑 구조가 비슷한데 더 싼 거 있어?"라고 물으면,
*   **Search**: 벡터 유사도(Cosine Similarity) + 가격 정렬(ASC)로 Ambeed 데이터를 추천.
