import os
import logging
from dotenv import load_dotenv
from app.core.supabase import supabase

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

load_dotenv(os.path.join(os.path.dirname(__file__), '.env'))

def check_columns():
    try:
        # Try to select the new columns
        # If they don't exist, this should fail
        response = supabase.table("commercial_reagents").select("id, ai_refined, target, properties, summary, source_name").limit(1).execute()
        logger.info("✅ Columns exist!")
        logger.info(f"Data sample: {response.data}")
    except Exception as e:
        logger.error(f"❌ Columns check failed: {e}")
        # Check error message for specific column missing
        if "column" in str(e) and "does not exist" in str(e):
             logger.info("⚠️ Columns are missing.")
        else:
             logger.info("⚠️ Error might be unrelated to columns.")

if __name__ == "__main__":
    check_columns()
