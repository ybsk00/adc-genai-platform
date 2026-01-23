"""
Final Report Query
"""
import os
import sys
from dotenv import load_dotenv
from supabase import create_client

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
load_dotenv()

supabase = create_client(
    os.getenv("SUPABASE_URL"),
    os.getenv("SUPABASE_SERVICE_KEY")
)

def run_report():
    print("📊 Running Final Report Query...")
    
    # Try to select extracted_target if it exists, otherwise just summary
    try:
        res = supabase.table("knowledge_base")\
            .select("title, summary, relevance_score")\
            .eq("source_type", "PubMed")\
            .gt("relevance_score", 0.7)\
            .order("created_at", desc=True)\
            .limit(5)\
            .execute()
            
        print(f"Found {len(res.data)} records.")
        for item in res.data:
            print("-" * 50)
            print(f"Title: {item['title'][:80]}...")
            print(f"Score: {item['relevance_score']}")
            print(f"Summary: {item['summary'][:200]}...")
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    run_report()
