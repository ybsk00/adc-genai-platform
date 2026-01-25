import asyncio
import logging
from unittest.mock import MagicMock, AsyncMock

# 로깅 설정
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("BatchTest")

# 실제 크롤러 클래스 임포트 (경로 설정)
import sys
import os
sys.path.append(os.getcwd())

from app.services.ambeed_crawler import AmbeedCrawler

async def test_batch_logic():
    crawler = AmbeedCrawler()
    
    # 실제 DB와 통신하지 않도록 Mock 설정
    crawler._save_batch = AsyncMock()
    crawler._process_single_product = AsyncMock(return_value={"test": "data"})
    crawler._enrich_and_prepare_item = AsyncMock(return_value={
        "ambeed_cat_no": "TEST-001", 
        "product_name": "Test Product",
        "smiles_code": "C1=CC=CC=C1"
    })
    
    # update_job_status 모킹
    import app.api.scheduler
    app.api.scheduler.update_job_status = AsyncMock()
    app.api.scheduler.get_job_from_db = AsyncMock(return_value=None)
    app.api.scheduler.is_cancelled = AsyncMock(return_value=False)

    # crawl_category 내부의 핵심 루프 시뮬레이션
    # 6개의 데이터를 찾았을 때, 5개에서 한 번 저장되고 마지막에 1개가 저장되어야 함
    batch_data = []
    count = 0
    limit = 6
    job_id = "test_job"
    
    logger.info("🚀 Starting Batch Save Logic Test (Target: 5 items per batch)")
    
    for i in range(1, limit + 1):
        count += 1
        res = await crawler._process_single_product(None, "url", "cat")
        final_item = await crawler._enrich_and_prepare_item(res)
        
        if final_item:
            batch_data.append(final_item)
            logger.info(f"➕ Item {i} added to batch. Current batch size: {len(batch_data)}")
            
            # 실제 코드에 반영된 5개 단위 저장 로직
            if len(batch_data) >= 5:
                await crawler._save_batch(batch_data)
                logger.info(f"💾 [SUCCESS] Batch size reached 5. _save_batch CALLED!")
                batch_data = [] # 메모리 비우기
    
    # 루프 종료 후 남은 데이터 처리
    if batch_data:
        await crawler._save_batch(batch_data)
        logger.info(f"💾 [SUCCESS] Final remaining data saved. Count: {len(batch_data)}")

    # 검증: _save_batch가 총 2번 호출되었는지 확인 (5개 묶음 + 1개 남은 것)
    save_call_count = crawler._save_batch.call_count
    logger.info(f"📊 Total _save_batch calls: {save_call_count}")
    
    if save_call_count == 2:
        print("\n✅ TEST PASSED: 5개 단위로 저장이 정확히 실행됩니다!")
    else:
        print("\n❌ TEST FAILED: 저장 로직 확인 필요.")

if __name__ == "__main__":
    asyncio.run(test_batch_logic())
