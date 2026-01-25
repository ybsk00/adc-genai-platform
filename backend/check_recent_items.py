import asyncio
from app.core.supabase import supabase
from datetime import datetime, timedelta

async def check_recent():
    # Check for items created in the last 5 minutes
    five_mins_ago = (datetime.utcnow() - timedelta(minutes=5)).isoformat()
    res = supabase.table("commercial_reagents").select("count", count="exact").gt("created_at", five_mins_ago).execute()
    print(f"New items in last 5 mins: {res.count}")

if __name__ == "__main__":
    asyncio.run(check_recent())
