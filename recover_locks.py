import os
import sys
from dotenv import load_dotenv

# 1. Load .env manually and inject into os.environ
env_path = os.path.join(os.getcwd(), "backend", ".env")
if os.path.exists(env_path):
    with open(env_path, "r") as f:
        for line in f:
            if "=" in line and not line.startswith("#"):
                key, value = line.strip().split("=", 1)
                os.environ[key] = value.strip("'").strip('"')
    print("Environment variables injected manually.")
else:
    print(f"Error: .env not found at {env_path}")
    sys.exit(1)

# 2. Add backend to sys.path
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.core.supabase import supabase

def recover_locks():
    print("--- Recovering Job Locks ---")
    
    # 3. Check if we have a real client
    if hasattr(supabase, "table") and supabase.__class__.__name__ == "MockSupabaseClient":
        print("Error: Still using MockSupabaseClient. Check your environment variables.")
        return

    # 4. Delete all job locks
    try:
        # Use a filter that matches all rows
        res = supabase.table("job_locks").delete().neq("job_type", "NONE_EXISTING_TYPE").execute()
        print(f"Deleted job locks: {len(res.data) if res.data else 0}")
    except Exception as e:
        print(f"Error deleting locks: {e}")
    
    # 5. Set running jobs to failed
    try:
        res = supabase.table("sync_jobs").update({
            "status": "failed",
            "errors": ["Force stopped due to lock recovery"],
            "completed_at": "2026-01-21T08:55:00Z"
        }).eq("status", "running").execute()
        print(f"Updated running jobs to failed: {len(res.data) if res.data else 0}")
    except Exception as e:
        print(f"Error updating jobs: {e}")
    
    print("\n--- Recovery Complete ---")

if __name__ == "__main__":
    recover_locks()
