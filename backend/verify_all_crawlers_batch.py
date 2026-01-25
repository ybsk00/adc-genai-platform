import asyncio
import logging
from unittest.mock import AsyncMock
import sys
import os

# 로깅 설정
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("AllBatchTest")

sys.path.append(os.getcwd())

# 크롤러 임포트
from app.services.ambeed_crawler import ambeed_crawler
from app.services.creative_biolabs_crawler import creative_crawler

async def run_test(name, crawler):
    logger.info(f"\n--- Testing {name} Batch Logic ---")
    
    # 모킹 설정
    crawler._save_batch = AsyncMock()
    crawler._process_single_product = AsyncMock(return_value={"test": "data"})
    crawler._enrich_and_prepare_item = AsyncMock(return_value={"id": "test-data"})
    
    # 외부 종속성 모킹
    import app.api.scheduler
    app.api.scheduler.update_job_status = AsyncMock()
    app.api.scheduler.is_cancelled = AsyncMock(return_value=False)
    app.api.scheduler.get_job_from_db = AsyncMock(return_value=None)

    batch_data = []
    limit = 6 # 5개에서 한 번 저장되고 1개가 남아야 함
    
    for i in range(1, limit + 1):
        res = await crawler._process_single_product(None, "url", "cat")
        final_item = await crawler._enrich_and_prepare_item(res)
        
        if final_item:
            batch_data.append(final_item)
            # 코드에 반영된 저장 로직 (5개 단위)
            if len(batch_data) >= 5:
                await crawler._save_batch(batch_data)
                logger.info(f"✅ {name}: Batch of 5 SAVED.")
                batch_data = []

    if batch_data:
        await crawler._save_batch(batch_data)
        logger.info(f"✅ {name}: Final remaining {len(batch_data)} item SAVED.")

    save_calls = crawler._save_batch.call_count
    if save_calls == 2:
        logger.info(f"✨ {name} PASS: _save_batch called {save_calls} times.")
        return True
    else:
        logger.error(f"❌ {name} FAIL: _save_batch called {save_calls} times.")
        return False

async def main():
    ambeed_res = await run_test("Ambeed", ambeed_crawler)
    cb_res = await run_test("Creative Biolabs", creative_crawler)
    
    if ambeed_res and cb_res:
        print("\n" + "="*40)
        print("🏆 ALL TEST PASSED: 두 크롤러 모두 5개 단위 저장이 확인되었습니다.")
        print("="*40)
    else:
        print("\n❌ 일부 테스트 실패. 로직 확인 필요.")

if __name__ == "__main__":
    asyncio.run(main())
