import os
import logging
from dotenv import load_dotenv

# Load .env explicitly
env_path = os.path.join(os.path.dirname(__file__), '.env')
print(f"Loading .env from: {env_path}")
load_dotenv(env_path)

from app.core.supabase import supabase

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_db():
    print(f"Supabase Client Type: {type(supabase)}")
    
    try:
        logger.info("🧪 Testing direct DB insert...")
        dummy_data = {
            "ambeed_cat_no": "TEST-DIRECT-002",
            "product_name": "Direct Insert Test 2",
            "source_name": "Creative Biolabs",
            "crawled_at": "2026-01-01T00:00:00"
        }
        res = supabase.table("commercial_reagents").upsert(dummy_data, on_conflict="ambeed_cat_no").execute()
        logger.info(f"✅ Direct insert successful! Data: {res.data}")
    except Exception as e:
        logger.error(f"❌ Direct insert failed: {e}")

if __name__ == "__main__":
    test_db()
