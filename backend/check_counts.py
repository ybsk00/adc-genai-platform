import os
import sys
from dotenv import load_dotenv

# Load .env BEFORE imports
env_path = os.path.join(os.getcwd(), "backend", ".env")
load_dotenv(env_path)

# Add backend to sys.path if running as script
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.core.supabase import supabase

def check_counts():
    print("Checking counts...")
    
    # Count False
    res_false = supabase.table("commercial_reagents").select("count", count="exact").eq("ai_refined", False).execute()
    count_false = res_false.count
    print(f"ai_refined = False: {count_false}")
    
    # Count NULL (Supabase client usually supports is_("column", "null"))
    try:
        res_null = supabase.table("commercial_reagents").select("count", count="exact").is_("ai_refined", "null").execute()
        count_null = res_null.count
        print(f"ai_refined = NULL: {count_null}")
    except Exception as e:
        print(f"Error checking NULLs: {e}")
        
    # Total
    res_total = supabase.table("commercial_reagents").select("count", count="exact").execute()
    print(f"Total rows: {res_total.count}")

    print("\n--- Golden Set Library ---")
    res_gs_false = supabase.table("golden_set_library").select("count", count="exact").eq("ai_refined", False).execute()
    print(f"ai_refined = False: {res_gs_false.count}")
    
    try:
        res_gs_null = supabase.table("golden_set_library").select("count", count="exact").is_("ai_refined", "null").execute()
        print(f"ai_refined = NULL: {res_gs_null.count}")
    except Exception as e:
        print(f"Error checking NULLs: {e}")
        
    res_gs_total = supabase.table("golden_set_library").select("count", count="exact").execute()
    print(f"Total rows: {res_gs_total.count}")

if __name__ == "__main__":
    check_counts()
