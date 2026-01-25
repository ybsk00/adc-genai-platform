import asyncio
import logging
import sys
import os

# backend 경로 추가
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.services.ambeed_crawler import ambeed_crawler
from app.core.supabase import supabase

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("TestSave")

async def test_save():
    job_id = "test_run_001"
    limit = 2
    start_page = 10  # 1페이지는 중복이 많으니 10페이지부터 테스트
    
    print(f"\n🚀 테스트 시작: Ambeed 수집 (Page {start_page}부터 {limit}개)")
    
    # 1. 실행 전 개수 확인
    res_before = supabase.table("commercial_reagents").select("count", count="exact").execute()
    count_before = res_before.count
    print(f"📊 수집 전 DB 레코드 수: {count_before}")
    
    # 2. 크롤러 실행
    # (API가 아닌 내부 메서드를 직접 호출하여 결과 확인)
    await ambeed_crawler.run(search_term="ADC Toxins", limit=limit, job_id=job_id, start_page=start_page)
    
    # 3. 실행 후 개수 확인
    res_after = supabase.table("commercial_reagents").select("count", count="exact").execute()
    count_after = res_after.count
    print(f"📊 수집 후 DB 레코드 수: {count_after}")
    
    # 4. 최신 저장 데이터 2개 출력
    latest = supabase.table("commercial_reagents")\
        .select("ambeed_cat_no, product_name, crawled_at")\
        .order("crawled_at", desc=True)\
        .limit(2)\
        .execute()
    
    print("\n✅ 최근 저장된 데이터 (DB 직접 조회):")
    for item in latest.data:
        print(f"- {item['ambeed_cat_no']}: {item['product_name']} (수집시간: {item['crawled_at']})")

if __name__ == "__main__":
    asyncio.run(test_save())
