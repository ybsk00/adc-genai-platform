import asyncio
import logging
import sys
from app.services.ambeed_crawler import ambeed_crawler

# Configure logging to stdout
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)

async def test_crawler():
    print("🚀 Starting Test Crawl (ADC Linker, Limit 5)...")
    
    # Use a category that likely has new items
    await ambeed_crawler.run(search_term='ADC Linker', limit=5, job_id='debug_test_1')
    
    print("✅ Test Crawl Finished")

if __name__ == "__main__":
    asyncio.run(test_crawler())
