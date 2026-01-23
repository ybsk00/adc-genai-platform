import os
import logging
import asyncio
from dotenv import load_dotenv

# 1. Load .env explicitly BEFORE any app imports
env_path = os.path.join(os.path.dirname(__file__), '.env')
print(f"Loading .env from: {env_path}")
load_dotenv(env_path)

# 2. Check env vars
print(f"SUPABASE_URL: {os.getenv('SUPABASE_URL')}")

# 3. Import app modules
from app.core.supabase import supabase
from app.services.creative_biolabs_crawler import creative_crawler

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def test_crawler():
    print(f"Supabase Client Type: {type(supabase)}")
    
    # Verify client is not Mock
    if "Mock" in str(type(supabase)):
        logger.error("❌ Using Mock Client! Aborting.")
        return

    logger.info("🧪 Testing Creative Biolabs Crawler...")
    
    category = "ADC Linker"
    url = creative_crawler.CATEGORIES[category]
    
    try:
        # Crawl 3 items
        count = await creative_crawler.crawl_category(category, url, limit=3)
        logger.info(f"✅ Crawled {count} items.")
    except Exception as e:
        logger.error(f"❌ Test failed: {e}")

if __name__ == "__main__":
    asyncio.run(test_crawler())
