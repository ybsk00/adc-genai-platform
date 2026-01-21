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

def check_failed_count():
    print("--- Checking Failed Records ---")
    
    try:
        if hasattr(supabase, "table") and supabase.__class__.__name__ == "MockSupabaseClient":
            print("Error: Still using MockSupabaseClient.")
            return

        # Check for 'LLM analysis failed' or any error
        res = supabase.table("golden_set_library")\
            .select("count", count="exact")\
            .neq("processing_error", "null")\
            .execute()
        
        total_errors = res.count
        print(f"Total records with processing_error: {total_errors}")
        
        if total_errors > 0:
            # Sample errors
            res = supabase.table("golden_set_library")\
                .select("processing_error")\
                .neq("processing_error", "null")\
                .limit(5)\
                .execute()
            print("Sample errors:")
            for item in res.data:
                print(f"- {item['processing_error']}")
                
    except Exception as e:
        print(f"Error checking DB: {e}")

if __name__ == "__main__":
    check_failed_count()
