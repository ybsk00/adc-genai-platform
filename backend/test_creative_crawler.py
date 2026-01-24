import asyncio
import logging
import sys
import os
from dotenv import load_dotenv

# Load .env explicitly
load_dotenv(os.path.join(os.path.dirname(__file__), '.env'))

# Add backend directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'app'))

from app.services.creative_biolabs_crawler import creative_crawler

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def test_crawler():
    logger.info("🧪 Starting Creative Biolabs Crawler Test...")
    
    # Test specific category: ADC Linker
    category = "ADC Linker"
    url = creative_crawler.CATEGORIES.get(category)
    
    if not url:
        logger.error(f"Category {category} not found!")
        return

    logger.info(f"Target URL: {url}")
    
    # Run crawler with limit 5 (small batch for testing)
    try:
        count = await creative_crawler.crawl_category(category, url, limit=5)
        logger.info(f"✅ Test Complete. Processed {count} items.")
    except Exception as e:
        logger.error(f"❌ Test Failed: {e}")

if __name__ == "__main__":
    asyncio.run(test_crawler())
