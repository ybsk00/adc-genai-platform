import asyncio
import logging
from app.services.ambeed_crawler import ambeed_crawler

# Configure logging
logging.basicConfig(level=logging.INFO)

async def test_crawler():
    print("🚀 Starting Test Crawl (Limit 2 per category)...")
    
    # Run for a specific category to be quick, or 'all' to test concurrency
    # Let's test 'ADC Toxins' specifically first, or 'all' with small limit
    await ambeed_crawler.run(search_term='ADC Toxins', limit=2, job_id='test_job_123')
    
    print("✅ Test Crawl Finished")

if __name__ == "__main__":
    asyncio.run(test_crawler())
