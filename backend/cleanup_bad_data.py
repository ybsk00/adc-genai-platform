"""
Data Cleanup Script
Deletes bad data from knowledge_base created after 10 AM today with score 0.
"""
import os
import sys
from dotenv import load_dotenv
from supabase import create_client

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

load_dotenv()

supabase = create_client(
    os.getenv("SUPABASE_URL"),
    os.getenv("SUPABASE_SERVICE_KEY")
)

def cleanup_bad_data():
    print("🧹 Starting cleanup of bad PubMed data...")
    
    # Check count before delete
    res = supabase.table("knowledge_base")\
        .select("id", count="exact")\
        .eq("source_type", "PubMed")\
        .eq("relevance_score", 0)\
        .gte("created_at", "2026-01-23T01:00:00")\
        .execute()
    
    count = res.count
    print(f"📋 Found {count} bad records to delete.")
    
    if count > 0:
        # Delete
        supabase.table("knowledge_base")\
            .delete()\
            .eq("source_type", "PubMed")\
            .eq("relevance_score", 0)\
            .gte("created_at", "2026-01-23T10:00:00")\
            .execute()
        print(f"✅ Deleted {count} records.")
    else:
        print("✨ No bad records found.")

if __name__ == "__main__":
    cleanup_bad_data()
