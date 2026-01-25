import asyncio
from app.core.supabase import supabase

async def check_count():
    res = supabase.table("commercial_reagents").select("count", count="exact").execute()
    print(f"Total items: {res.count}")

if __name__ == "__main__":
    asyncio.run(check_count())
