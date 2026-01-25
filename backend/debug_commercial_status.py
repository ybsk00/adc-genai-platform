import asyncio
import os
from app.core.supabase import supabase

async def check_status():
    print("Checking commercial_reagents status...")
    
    # 1. Total count
    res = supabase.table("commercial_reagents").select("count", count="exact").execute()
    total = res.count
    print(f"Total records: {total}")
    
    # 2. Refined count
    res = supabase.table("commercial_reagents").select("count", count="exact").eq("ai_refined", True).execute()
    refined = res.count
    print(f"Refined records: {refined}")
    
    # 3. Pending count
    res = supabase.table("commercial_reagents").select("count", count="exact").eq("ai_refined", False).execute()
    pending = res.count
    print(f"Pending records: {pending}")

if __name__ == "__main__":
    asyncio.run(check_status())
