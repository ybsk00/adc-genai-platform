import asyncio
import os
from app.core.supabase import supabase

async def reset_status():
    print("Resetting commercial_reagents status to pending...")
    
    # Update all records where ai_refined is True to False
    # We can just update all records since the user wants to restart
    res = supabase.table("commercial_reagents").update({"ai_refined": False}).neq("id", "00000000-0000-0000-0000-000000000000").execute()
    
    # Check how many were updated (Supabase update returns the data)
    if res.data:
        print(f"Successfully reset {len(res.data)} records.")
    else:
        print("No records updated or error occurred.")

    # Verify
    res = supabase.table("commercial_reagents").select("count", count="exact").eq("ai_refined", False).execute()
    print(f"Pending records after reset: {res.count}")

if __name__ == "__main__":
    asyncio.run(reset_status())
