"""
Creative Biolabs 크롤러 테스트 스크립트
"""
import aiohttp
import asyncio
from bs4 import BeautifulSoup

async def test_creative_biolabs():
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'
    }
    
    # 테스트 1: HER2 검색
    print("=== Test 1: HER2 Search ===")
    async with aiohttp.ClientSession(headers=headers) as session:
        url = "https://www.creative-biolabs.com/adc/search.aspx?q=HER2"
        async with session.get(url) as res:
            print(f"Status: {res.status}")
            if res.status == 200:
                html = await res.text()
                soup = BeautifulSoup(html, "html.parser")
                print(f"Title: {soup.title.string if soup.title else 'No title'}")
                
                # 제품 링크 찾기
                product_links = soup.select("a[href*='/adc/']")
                print(f"Found {len(product_links)} links containing '/adc/'")
                for link in product_links[:5]:
                    print(f"  - {link.get('href')}")
            else:
                print(f"Error: {await res.text()[:500]}")

    # 테스트 2: 카테고리 페이지
    print("\n=== Test 2: Category Page ===")
    async with aiohttp.ClientSession(headers=headers) as session:
        url = "https://www.creative-biolabs.com/adc/products.html"
        async with session.get(url) as res:
            print(f"Status: {res.status}")
            if res.status == 200:
                html = await res.text()
                soup = BeautifulSoup(html, "html.parser")
                print(f"Title: {soup.title.string if soup.title else 'No title'}")
                
                # 제품 링크 찾기
                product_links = soup.select("a[href*='/adc/p/'], a[href*='/adc/product/']")
                print(f"Found {len(product_links)} product links")
                for link in product_links[:5]:
                    print(f"  - {link.get('href')}")
            else:
                print(f"Error: {res.status}")

    # 테스트 3: 모든 링크 구조 분석
    print("\n=== Test 3: All ADC Links Analysis ===")
    async with aiohttp.ClientSession(headers=headers) as session:
        url = "https://www.creative-biolabs.com/adc/"
        async with session.get(url) as res:
            print(f"Status: {res.status}")
            if res.status == 200:
                html = await res.text()
                soup = BeautifulSoup(html, "html.parser")
                
                all_links = soup.select("a[href*='/adc/']")
                unique_patterns = set()
                for link in all_links:
                    href = link.get('href', '')
                    # 패턴 추출
                    parts = href.split('/')
                    if len(parts) > 2:
                        pattern = '/'.join(parts[:3])
                        unique_patterns.add(pattern)
                
                print(f"Unique link patterns: {len(unique_patterns)}")
                for p in sorted(unique_patterns)[:10]:
                    print(f"  - {p}")

asyncio.run(test_creative_biolabs())
