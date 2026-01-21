import os
import sys
from dotenv import load_dotenv
import traceback

# Load env
load_dotenv("backend/.env")

# Add backend to path
sys.path.append("backend")

from app.core.supabase import supabase
from app.core.config import settings
from app.services.ai_refiner import ai_refiner
import asyncio
import json

async def main():
    print("=== [1] Clinical Trials (Gemini 2.0 Flash - Direct SDK) Test ===")
    from app.services.ai_refiner import ai_refiner
    
    # Get a failed record
    res = supabase.table("golden_set_library")\
        .select("*")\
        .not_.is_("processing_error", "null")\
        .limit(1)\
        .execute()
    
    if not res.data:
        res = supabase.table("golden_set_library").select("*").eq("ai_refined", False).limit(1).execute()
    
    if res.data:
        record = res.data[0]
        print(f"Record ID: {record['id']}")
        print(f"Record Name: {record.get('name')}")
        result = await ai_refiner.refine_single_record(record)
        print(f"RESULT:\n{json.dumps(result, indent=2)}")
    else:
        print("No records found for Clinical Trials test.")
    
    print("\n=== [2] PubMed (Gemini 2.0 Flash - Direct SDK) Test ===")
    from app.services.knowledge_refiner import knowledge_refiner
    
    # Get a PubMed record
    res = supabase.table("knowledge_base")\
        .select("*")\
        .limit(1)\
        .execute()
    
    if res.data:
        item = res.data[0]
        print(f"Item ID: {item['id']}")
        result = await knowledge_refiner.analyze_abstract(item["content"])
        print(f"RESULT:\n{json.dumps(result, indent=2)}")
    else:
        print("No PubMed records found.")

    print("\n=== [3] Cost Tracking Verification ===")
    res = supabase.table("llm_usage_logs").select("*").order("created_at", desc=True).limit(5).execute()
    print(f"Recent Logs:\n{json.dumps(res.data, indent=2)}")

if __name__ == "__main__":
    asyncio.run(main())
