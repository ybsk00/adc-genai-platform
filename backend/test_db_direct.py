import asyncio
import logging
import sys
import os

# 프로젝트 루트 경로 추가
sys.path.append(os.getcwd())

from app.core.supabase import supabase
from datetime import datetime

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DB_Test")

async def test_write():
    test_id = f"TEST-LOG-{datetime.now().strftime('%H%M%S')}"
    logger.info(f"🚀 [팩트체크] DB 쓰기 테스트 시작 (ID: {test_id})")
    
    test_data = {
        "ambeed_cat_no": test_id,
        "product_name": "DB WRITE TEST RECORD",
        "category": "TEST",
        "source_name": "TEST_RUNNER",
        "crawled_at": datetime.utcnow().isoformat()
    }
    
    try:
        # 1. UPSERT 시도
        logger.info(f"📡 Supabase에 데이터 전송 중... (URL: {os.getenv('SUPABASE_URL')})")
        res = supabase.table("commercial_reagents").upsert(test_data, on_conflict="ambeed_cat_no").execute()
        
        if res.data:
            logger.info("✅ [성공] DB에 데이터가 정상적으로 저장되었습니다!")
            print(f"\n[RESULT] DB_WRITE_SUCCESS: {res.data[0]['id']}")
        else:
            logger.error("❌ [실패] 응답 데이터가 없습니다. 권한 설정을 확인하세요.")
            
    except Exception as e:
        logger.error(f"🔥 [치명적 에러] DB 접속 중 오류 발생: {e}")
        print(f"\n[RESULT] DB_WRITE_FAILED: {str(e)}")

if __name__ == "__main__":
    asyncio.run(test_write())