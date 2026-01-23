import asyncio
import aiohttp
import sys

API_BASE_URL = "http://localhost:8000"

async def reproduce_error():
    async with aiohttp.ClientSession() as session:
        print("🚀 Triggering Creative Biolabs Crawler (Category=All, Limit=1000)...")
        # The user reported: run?category=all&limit=1000
        url = f"{API_BASE_URL}/api/scheduler/crawler/creative-biolabs/run?category=all&limit=1000"
        
        try:
            async with session.post(url) as resp:
                print(f"Status Code: {resp.status}")
                text = await resp.text()
                print(f"Response: {text}")
                
                if resp.status == 500:
                    print("❌ Reproduced 500 Error!")
                else:
                    print("✅ Request succeeded (or at least didn't return 500 immediately).")
        except Exception as e:
            print(f"❌ Request failed: {e}")

if __name__ == "__main__":
    asyncio.run(reproduce_error())
