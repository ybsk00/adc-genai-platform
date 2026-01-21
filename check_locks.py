import os
import json
from dotenv import load_dotenv

# Load .env before any app imports
env_path = os.path.join(os.getcwd(), "backend", ".env")
load_dotenv(env_path)

# Add backend to sys.path
import sys
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.core.supabase import supabase

def check_locks():
    print(f"--- Environment Check ---")
    print(f"SUPABASE_URL: {os.getenv('SUPABASE_URL')[:15]}...")
    
    print("\n--- Current Job Locks ---")
    res = supabase.table("job_locks").select("*").execute()
    print(json.dumps(res.data, indent=2))
    
    print("\n--- Recent Sync Jobs (Running) ---")
    res = supabase.table("sync_jobs").select("*").eq("status", "running").order("started_at", desc=True).limit(5).execute()
    print(json.dumps(res.data, indent=2))

if __name__ == "__main__":
    check_locks()
