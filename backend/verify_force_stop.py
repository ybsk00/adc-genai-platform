import asyncio
import aiohttp
import time
from app.core.config import settings
from app.core.supabase import supabase

API_BASE_URL = "http://localhost:8000"

async def test_force_stop():
    async with aiohttp.ClientSession() as session:
        # 1. Start a long-running job (e.g., Clinical Trials Import)
        print("🚀 Starting Clinical Trials Import...")
        async with session.post(f"{API_BASE_URL}/api/scheduler/bulk/import?max_studies=1000&mode=full") as resp:
            if resp.status != 200:
                print(f"❌ Failed to start job: {await resp.text()}")
                return
            data = await resp.json()
            job_id = data["job_id"]
            print(f"✅ Job started: {job_id}")

        # 2. Wait a bit to let it run
        print("⏳ Waiting for 5 seconds...")
        await asyncio.sleep(5)

        # 3. Check status (should be running)
        async with session.get(f"{API_BASE_URL}/api/scheduler/sync/{job_id}") as resp:
            status_data = await resp.json()
            print(f"📊 Current Status: {status_data['status']}")

        # 4. Trigger Force Stop All Workers
        print("🛑 Triggering Force Stop All Workers...")
        async with session.post(f"{API_BASE_URL}/api/scheduler/workers/reset") as resp:
            if resp.status == 200:
                print("✅ Force stop request sent.")
            else:
                print(f"❌ Force stop failed: {await resp.text()}")

        # 5. Check status again (should be stopped)
        print("⏳ Waiting for 2 seconds...")
        await asyncio.sleep(2)
        
        async with session.get(f"{API_BASE_URL}/api/scheduler/sync/{job_id}") as resp:
            final_status = await resp.json()
            print(f"📊 Final Status: {final_status['status']}")
            print(f"📝 Errors: {final_status.get('errors')}")
            
            # Check if cancel_requested was set (we can't check this easily via API unless we query DB directly, but status 'stopped' implies it worked if the worker respected it)
            # However, reset_all_workers sets status to stopped immediately.
            # Ideally we should see logs in the backend console saying "Import cancelled by user"

if __name__ == "__main__":
    asyncio.run(test_force_stop())
