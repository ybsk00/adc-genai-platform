import asyncio
import os
import sys
from datetime import datetime

# Add backend directory to sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), 'app'))

from app.core.supabase import supabase

async def debug_db():
    with open("debug_result.txt", "w", encoding="utf-8") as f:
        f.write("--- Debugging DB Status ---\n")
        
        # 1. Check sync_jobs
        f.write("\n[Sync Jobs (Last 5)]\n")
        try:
            res = supabase.table("sync_jobs").select("*").order("started_at", desc=True).limit(5).execute()
            for job in res.data:
                f.write(f"ID: {job.get('id')}, Source: {job.get('source')}, Status: {job.get('status')}, Started: {job.get('started_at')}, Errors: {job.get('errors')}\n")
        except Exception as e:
            f.write(f"Error checking sync_jobs: {e}\n")

        # 2. Check commercial_reagents count and sample
        f.write("\n[Commercial Reagents]\n")
        try:
            count_res = supabase.table("commercial_reagents").select("count", count="exact").execute()
            f.write(f"Total Count: {count_res.count}\n")
            
            # Check for recent entries
            res = supabase.table("commercial_reagents").select("*").order("crawled_at", desc=True).limit(3).execute()
            for item in res.data:
                f.write(f"ID: {item.get('id')}, Name: {item.get('product_name')}, AI Refined: {item.get('ai_refined')}, Target: {item.get('target')}\n")
                
        except Exception as e:
            f.write(f"Error checking commercial_reagents: {e}\n")

        # 3. Check LLM Usage
        f.write("\n[LLM Usage Today]\n")
        try:
            today = datetime.now().date().isoformat()
            res = supabase.table("llm_usage_logs").select("cost_usd").gte("created_at", f"{today}T00:00:00").execute()
            total_cost = sum(item['cost_usd'] for item in res.data)
            f.write(f"Total Cost Today (USD): {total_cost}\n")
        except Exception as e:
            f.write(f"Error checking LLM usage: {e}\n")

if __name__ == "__main__":
    asyncio.run(debug_db())
