import asyncio
import os
import sys
import logging
from dotenv import load_dotenv

# backend 경로 추가
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.append(current_dir)

# 환경변수 로드
load_dotenv(os.path.join(current_dir, ".env"))

from app.services.ambeed_crawler import ambeed_crawler
from app.core.supabase import supabase

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("BulkScraper")

async def run_bulk():
    # 사장님 지시사항 반영
    # 1. 2개씩 하지 말고 크게 묶어서 (batch_size=50)
    # 2. 최대한 많이 (limit=5000)
    # 3. 51페이지부터 (과거 데이터 구역)
    
    start_page = 51
    limit = 5000
    batch_size = 50
    job_id = "bulk_manual_run_001"
    
    print(f"\n🔥 [대량 수집 시작] 시작페이지: {start_page}, 목표량: {limit}, 저장단위: {batch_size}")
    
    try:
        # DB 연결 확인
        res = supabase.table("commercial_reagents").select("count", count="exact").limit(1).execute()
        print(f"✅ DB 연결 확인됨. 현재 데이터 수: {res.count}")
        
        # 크롤러 실행
        await ambeed_crawler.run(
            search_term="all", 
            limit=limit, 
            job_id=job_id, 
            start_page=start_page,
            batch_size=batch_size
        )
        
    except Exception as e:
        logger.error(f"❌ 대량 수집 중 에러 발생: {e}", exc_info=True)

if __name__ == "__main__":
    asyncio.run(run_bulk())
