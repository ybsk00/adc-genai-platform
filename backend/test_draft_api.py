import asyncio
import os
import sys
from dotenv import load_dotenv

# Load .env
env_path = os.path.join(os.getcwd(), "backend", ".env")
load_dotenv(env_path)
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.core.supabase import supabase

async def test_draft_api():
    print("Testing get_golden_set_drafts logic...")
    try:
        limit = 10
        offset = 0
        search = ""
        source = ""
        
        count_query = supabase.table("golden_set_library").select("*", count="exact").eq("status", "draft")
        count_res = count_query.execute()
        print(f"Count: {count_res.count}")
        
        data_query = supabase.table("golden_set_library").select("*").eq("status", "draft")
        response = data_query.order("created_at", desc=True).range(offset, offset + limit - 1).execute()
        
        print(f"Fetched {len(response.data)} items.")
        if response.data:
            print("Sample item:", response.data[0]['name'])
            
            # Simulate mapping
            drafts = []
            for item in response.data:
                try:
                    drafts.append({
                        "id": item.get("id"),
                        "drug_name": item.get("name") or "Unknown",
                        "target": item.get("properties", {}).get("target") or "Unknown",
                        # ... other fields
                    })
                except Exception as e:
                    print(f"Error mapping item {item.get('id')}: {e}")
            print(f"Mapped {len(drafts)} items successfully.")
            
    except Exception as e:
        print(f"API Logic Error: {e}")

if __name__ == "__main__":
    asyncio.run(test_draft_api())
