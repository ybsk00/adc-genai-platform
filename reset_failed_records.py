import os
import sys
from dotenv import load_dotenv

# 1. Load .env manually
env_path = os.path.join(os.getcwd(), "backend", ".env")
if os.path.exists(env_path):
    with open(env_path, "r") as f:
        for line in f:
            if "=" in line and not line.startswith("#"):
                key, value = line.strip().split("=", 1)
                os.environ[key] = value.strip("'").strip('"')
else:
    print("Error: .env not found")
    sys.exit(1)

# 2. Add backend to sys.path
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.core.supabase import supabase

def reset_failed_records():
    print("--- Resetting Failed Records ---")
    
    try:
        # Update records where processing_error is not null
        # Set ai_refined = false, processing_error = null
        res = supabase.table("golden_set_library")\
            .update({
                "ai_refined": False,
                "processing_error": None,
                "outcome_type": "Unknown"
            })\
            .neq("processing_error", "null")\
            .execute()
        
        print(f"Reset complete. Affected rows: {len(res.data) if res.data else 0}")
                
    except Exception as e:
        print(f"Error resetting DB: {e}")

if __name__ == "__main__":
    reset_failed_records()
