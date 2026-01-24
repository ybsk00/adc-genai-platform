import os
import sys
from datetime import datetime

# Add backend directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'app'))

from app.core.supabase import supabase

def test_insert():
    print("🧪 Testing DB Insert...")
    
    # 1. Try inserting without source_name
    data_v1 = {
        "ambeed_cat_no": "TEST-001",
        "product_name": "Test Product 1",
        "crawled_at": datetime.utcnow().isoformat()
    }
    
    try:
        res = supabase.table("commercial_reagents").upsert(data_v1, on_conflict="ambeed_cat_no").execute()
        print("✅ Insert V1 (no source_name) Success!")
    except Exception as e:
        print(f"❌ Insert V1 Failed: {e}")

    # 2. Try inserting with source_name
    data_v2 = {
        "ambeed_cat_no": "TEST-002",
        "product_name": "Test Product 2",
        "source_name": "Creative Biolabs",
        "crawled_at": datetime.utcnow().isoformat()
    }
    
    try:
        res = supabase.table("commercial_reagents").upsert(data_v2, on_conflict="ambeed_cat_no").execute()
        print("✅ Insert V2 (with source_name) Success!")
    except Exception as e:
        print(f"❌ Insert V2 Failed: {e}")

if __name__ == "__main__":
    with open("db_test_result.txt", "w", encoding="utf-8") as f:
        sys.stdout = f
        test_insert()
