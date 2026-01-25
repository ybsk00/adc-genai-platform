import asyncio
import logging
import sys
import os

# Ensure backend directory is in python path
sys.path.append(os.path.join(os.path.dirname(__file__)))

from app.services.ambeed_crawler import ambeed_crawler

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger("LocalRunner")

async def main():
    logger.info("🚀 Starting Local Ambeed Crawler Test (1 Item)...")
    
    # Check if we can import scheduler or if we need to mock it
    try:
        from app.api.scheduler import update_job_status
        logger.info("✅ app.api.scheduler found.")
    except ImportError:
        logger.warning("⚠️ app.api.scheduler not found. Mocking it to prevent crash.")
        from unittest.mock import MagicMock
        sys.modules["app.api.scheduler"] = MagicMock()
        sys.modules["app.api.scheduler"].update_job_status = asyncio.coroutine(lambda *args, **kwargs: None)
        sys.modules["app.api.scheduler"].is_cancelled = asyncio.coroutine(lambda *args, **kwargs: False)
        sys.modules["app.api.scheduler"].get_job_from_db = asyncio.coroutine(lambda *args, **kwargs: {})

    try:
        # Run for 1 item
        await ambeed_crawler.run(
            search_term="ADC Toxins", 
            limit=1, 
            job_id="local_test_run_001",
            start_page=1,
            batch_size=1
        )
        logger.info("✅ Test Run Completed Successfully.")
    except Exception as e:
        logger.error(f"❌ Test Run Failed: {e}", exc_info=True)

if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        logger.info("🛑 Stopped by user.")
