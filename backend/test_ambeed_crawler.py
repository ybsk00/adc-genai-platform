import asyncio
import logging
from app.services.ambeed_crawler import ambeed_crawler

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("test_crawler_error.log", encoding='utf-8'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

async def test_crawler():
    logger.info("Starting Ambeed Crawler Test...")
    
    # Test parameters
    search_term = "ADC Toxins" # Should match a key in CATEGORIES
    limit = 1 # Reduce limit for quick debug
    job_id = "test_job_123"
    
    try:
        # Run crawler
        await ambeed_crawler.run(search_term, limit, job_id)
        logger.info("Test completed successfully.")
    except Exception as e:
        logger.error(f"Test failed: {e}", exc_info=True)

if __name__ == "__main__":
    asyncio.run(test_crawler())
