import asyncio
import os
import sys
from pprint import pprint
from dotenv import load_dotenv

# Load .env explicitly
load_dotenv(os.path.join(os.path.dirname(__file__), '.env'))

# Add backend directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'app'))

from app.core.supabase import supabase

def verify_db():
    print("🔍 Verifying DB Results...")
    try:
        res = supabase.table("commercial_reagents")\
            .select("product_name, cas_number, smiles_code, target, properties, summary, price_data, product_url, crawled_at")\
            .order("crawled_at", desc=True)\
            .limit(5)\
            .execute()
            
        if not res.data:
            print("❌ No data found for Creative Biolabs!")
            return

        print(f"✅ Found {len(res.data)} records:")
        for item in res.data:
            print("-" * 50)
            pprint(item)
            
    except Exception as e:
        print(f"❌ DB Verification Failed: {e}")

if __name__ == "__main__":
    with open("db_verify_result.txt", "w", encoding="utf-8") as f:
        sys.stdout = f
        verify_db()
