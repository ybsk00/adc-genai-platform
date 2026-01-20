import os
from dotenv import load_dotenv
from supabase import create_client, Client

load_dotenv()

url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_KEY")

def check_columns():
    print("=== Checking Table Columns ===")
    supabase: Client = create_client(url, key)
    
    # Try to insert a dummy record with 'raw_data' to trigger the specific error
    # or just try to select it.
    # Unfortunately, Supabase-py doesn't have a direct 'get_columns' method easily accessible 
    # without using rpc or direct SQL if enabled.
    # But we can infer it from the error message we already got.
    
    print("Attempting to select 'raw_data' column...")
    try:
        supabase.table("golden_set_library").select("raw_data").limit(1).execute()
        print("SUCCESS: 'raw_data' column exists.")
    except Exception as e:
        print(f"FAILURE: Could not select 'raw_data'.")
        print(f"Error Details: {e}")
        print("\nCONCLUSION: The 'raw_data' column is MISSING from the table.")

if __name__ == "__main__":
    check_columns()
