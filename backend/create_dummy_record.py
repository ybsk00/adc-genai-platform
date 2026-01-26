from app.core.supabase import supabase
import asyncio

async def create_dummy():
    print("Creating dummy record...")
    res = supabase.table("golden_set_library").insert({
        "name": "Pazopanib Test Record Full",
        "properties": {}
    }).execute()
    print("Created:", res.data)

if __name__ == "__main__":
    asyncio.run(create_dummy())
