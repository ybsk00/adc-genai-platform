import asyncio
import os
import sys
from datetime import datetime

# Add backend directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), 'app'))

from app.core.supabase import supabase

async def debug_db():
    with open("debug_result_v2.txt", "w", encoding="utf-8") as f:
        f.write("--- Debugging DB Status V2 ---\n")
        
        # 1. Check specific job status
        job_id = "crawl_ambeed_6179af2b"
        f.write(f"\n[Checking Job: {job_id}]\n")
        try:
            res = supabase.table("sync_jobs").select("*").eq("id", job_id).execute()
            if res.data:
                job = res.data[0]
                f.write(f"Status: {job.get('status')}\n")
                f.write(f"Started: {job.get('started_at')}\n")
                f.write(f"Records Drafted: {job.get('records_drafted')}\n")
                f.write(f"Errors: {job.get('errors')}\n")
            else:
                f.write("Job not found.\n")
        except Exception as e:
            f.write(f"Error checking job: {e}\n")

        # 2. Check recent commercial_reagents (last 10)
        f.write("\n[Recent Commercial Reagents (Last 10)]\n")
        try:
            res = supabase.table("commercial_reagents").select("*").order("crawled_at", desc=True).limit(10).execute()
            for item in res.data:
                f.write(f"ID: {item.get('id')}\n")
                f.write(f"  Name: {item.get('product_name')}\n")
                f.write(f"  Crawled At: {item.get('crawled_at')}\n")
                f.write(f"  AI Refined: {item.get('ai_refined')}\n")
                f.write(f"  Target: {item.get('target')}\n")
                f.write("-" * 20 + "\n")
                
        except Exception as e:
            f.write(f"Error checking commercial_reagents: {e}\n")

if __name__ == "__main__":
    asyncio.run(debug_db())
