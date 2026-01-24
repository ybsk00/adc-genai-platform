import asyncio
import logging
import os
from dotenv import load_dotenv
from app.core.supabase import supabase
from app.services.ai_refiner import ai_refiner

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

load_dotenv(os.path.join(os.path.dirname(__file__), '.env'))

# Debug Env
print(f"CWD: {os.getcwd()}")
print(f"Env file path: {os.path.join(os.path.dirname(__file__), '.env')}")
print(f"SUPABASE_URL: {os.getenv('SUPABASE_URL')}")
print(f"SUPABASE_SERVICE_KEY: {os.getenv('SUPABASE_SERVICE_KEY')}")

async def test_upsert_and_trigger():
    logger.info("🚀 Testing Upsert and AI Trigger...")
    
    # 1. Insert a dummy record to commercial_reagents
    dummy_data = {
        "ambeed_cat_no": "TEST-DEBUG-001",
        "product_name": "Debug Test Reagent (MMAE)",
        "source_name": "Ambeed",
        "product_url": "https://example.com/debug",
        "category": "ADC Toxins",
        "ai_refined": False
    }
    
    try:
        # Mimic the crawler code
        # Note: supabase-py v2 upsert returns a PostgrestBuilder, need .select().execute()
        res = supabase.table("commercial_reagents").upsert(dummy_data, on_conflict="ambeed_cat_no").select().execute()
        
        logger.info(f"Response: {res}")
        if res.data:
            record = res.data[0]
            logger.info(f"✅ Upsert returned data: ID={record.get('id')}")
            
            # 2. Manually trigger refinement
            logger.info("🚀 Triggering Refinement manually...")
            analysis = await ai_refiner.refine_single_record(record)
            
            if analysis and "error" not in analysis:
                logger.info(f"✨ Refinement Success: {analysis}")
                
                # Update DB
                update_data = {
                    "target": analysis.get("target"),
                    "ai_refined": True,
                    "properties": {"ai_analysis": analysis}
                }
                supabase.table("commercial_reagents").update(update_data).eq("id", record["id"]).execute()
                logger.info("✅ DB Updated with AI result.")
            else:
                logger.error(f"❌ Refinement Failed: {analysis}")
                
        else:
            logger.error("❌ Upsert did NOT return data. Check .select() chaining.")
            
    except Exception as e:
        logger.error(f"❌ Exception during test: {e}")

if __name__ == "__main__":
    asyncio.run(test_upsert_and_trigger())
