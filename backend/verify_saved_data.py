import os
import logging
from dotenv import load_dotenv
from supabase import create_client

load_dotenv()
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("VerifyDB")

def verify():
    url = os.getenv("SUPABASE_URL")
    key = os.getenv("SUPABASE_SERVICE_KEY")
    
    if not url or not key:
        logger.error("Missing Supabase credentials.")
        return

    client = create_client(url, key)
    
    # Target Cat No from the previous log
    target_cat = "7-aminomethyl-10-methyl-11-fluoro-camptothecin"
    
    logger.info(f"🔍 Checking DB for: {target_cat}")
    
    res = client.table("commercial_reagents").select("*").eq("ambeed_cat_no", target_cat).execute()
    
    if res.data:
        logger.info("✅ Record FOUND in Database!")
        logger.info(f"   ID: {res.data[0]['id']}")
        logger.info(f"   Name: {res.data[0]['product_name']}")
        logger.info(f"   Created At: {res.data[0].get('crawled_at')}")
    else:
        logger.error("❌ Record NOT FOUND. The previous save might have been rolled back or blocked.")

if __name__ == "__main__":
    verify()
