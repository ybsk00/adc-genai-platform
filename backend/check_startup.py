import asyncio
import os
from app.core.supabase import supabase

async def check_and_fix_status():
    print("--- System Status Check ---")
    # 1. AI Refiner Status 확인
    res = supabase.table("system_config").select("*").eq("key", "AI_REFINER_STATUS").execute()
    if res.data:
        status = res.data[0]['value']
        print(f"Current AI_REFINER_STATUS: {status}")
        if status == "PAUSED":
            print("Action: Changing status to ACTIVE...")
            supabase.table("system_config").update({"value": "ACTIVE"}).eq("key", "AI_REFINER_STATUS").execute()
            print("Status updated to ACTIVE.")
    else:
        print("AI_REFINER_STATUS not found. Creating as ACTIVE...")
        supabase.table("system_config").insert({"key": "AI_REFINER_STATUS", "value": "ACTIVE"}).execute()

    # 2. Daily Cost 확인
    from app.services.cost_tracker import cost_tracker
    usage = await cost_tracker.get_usage_summary()
    print(f"Daily Usage: ${usage['daily_usage_usd']:.4f} / ${usage['daily_limit_usd']:.2f}")
    
    # 3. Pending Records 확인
    res = supabase.table("commercial_reagents").select("count", count="exact").eq("ai_refined", False).execute()
    print(f"Pending Commercial Reagents: {res.count}")

if __name__ == "__main__":
    asyncio.run(check_and_fix_status())