import os
import sys
import asyncio

# Load .env manually
env_path = os.path.join(os.getcwd(), "backend", ".env")
if os.path.exists(env_path):
    with open(env_path, "r", encoding="utf-8") as f:
        for line in f:
            if "=" in line and not line.strip().startswith("#"):
                key, value = line.strip().split("=", 1)
                # Handle quotes if present
                value = value.strip().strip("'").strip('"')
                os.environ[key] = value
else:
    print("Warning: .env not found at", env_path)

# Setup path to include 'backend' so we can import app modules
# Assuming we are running this from the project root
current_dir = os.getcwd()
backend_dir = os.path.join(current_dir, "backend")
if backend_dir not in sys.path:
    sys.path.append(backend_dir)

from app.core.supabase import supabase

async def reset_queue():
    job_id = "crawl_creative_c23157d0"
    print(f"Attempting to reset job: {job_id}")
    
    try:
        # Check current status
        # Note: We use .execute() directly as per the fix we just applied to the crawler, 
        # assuming standard usage is okay here.
        res = supabase.table("sync_jobs").select("*").eq("id", job_id).execute()
        
        if res.data:
            current_status = res.data[0].get('status')
            print(f"Current status: {current_status}")
            
            # Reset to pending
            update_res = supabase.table("sync_jobs").update({
                "status": "pending",
                "errors": None, 
                "completed_at": None
            }).eq("id", job_id).execute()
            
            print(f"✅ Job {job_id} reset to 'pending'. It should be picked up by the scheduler.")
        else:
            print(f"⚠️ Job {job_id} not found in sync_jobs table.")
            
    except Exception as e:
        print(f"❌ Error resetting queue: {e}")

if __name__ == "__main__":
    asyncio.run(reset_queue())
