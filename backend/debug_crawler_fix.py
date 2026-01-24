import os
import sys
import asyncio
import logging
from datetime import datetime

# 1. Load .env manually to ensure DB connection works
env_path = os.path.join(os.getcwd(), "backend", ".env")
if os.path.exists(env_path):
    with open(env_path, "r", encoding="utf-8") as f:
        for line in f:
            if "=" in line and not line.strip().startswith("#"):
                key, value = line.strip().split("=", 1)
                value = value.strip().strip("'").strip('"')
                os.environ[key] = value
else:
    print("⚠️ Warning: .env not found")

# 2. Add backend to sys.path
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.core.supabase import supabase

# Mock Logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DEBUG_TEST")

async def test_db_save_fix():
    print("🧪 Starting Debug Test for Crawler DB Fix...")
    
    # 3. Create dummy data (mimicking CreativeBiolabsCrawler output)
    test_data = {
        "ambeed_cat_no": "DEBUG-TEST-FIX-001",  # Unique ID for test
        "cas_number": "00-00-0",
        "product_name": "Debug Test Product (Val-Cit-PAB-Test)",
        "product_url": "https://example.com/debug-test",
        "category": "ADC Linker",
        "smiles_code": "C=C",
        "target": "Test Target",
        "properties": {"test": "true", "purity": ">99%"},
        "summary": "This is a debug test record to verify the DB save fix.",
        "price_data": [{"size": "1mg", "price": "$100", "availability": "In Stock"}],
        "source_name": "Creative Biolabs",
        "crawled_at": datetime.utcnow().isoformat()
    }

    print(f"   📝 Preparing to upsert record: {test_data['ambeed_cat_no']}")

    # 4. Try to save using the FIXED syntax
    # Old (Broken): .upsert(...).select().execute()
    # New (Fixed):  .upsert(...).execute()
    
    try:
        # Intentionally using the exact same line as the fix in creative_biolabs_crawler.py
        res = supabase.table("commercial_reagents").upsert(test_data, on_conflict="ambeed_cat_no").execute()
        
        print("\n✅ DB Save Successful!")
        print(f"   Response Data: {res.data}")
        
        if res.data:
            print(f"   Created/Updated ID: {res.data[0].get('id')}")
            print("   ✨ The '.select()' error is resolved. Data is safely landing in Supabase.")
            
    except Exception as e:
        print("\n❌ DB Save Failed!")
        print(f"   Error: {e}")
        print("   ⚠️ The fix might not be working or there is another issue.")

if __name__ == "__main__":
    asyncio.run(test_db_save_fix())
