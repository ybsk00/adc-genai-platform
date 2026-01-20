import os
from dotenv import load_dotenv
from supabase import create_client, Client

load_dotenv()

url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_KEY")

def check_table():
    print("=== Checking 'data_sync_logs' Table ===")
    supabase: Client = create_client(url, key)
    
    try:
        # Try to select from the table
        res = supabase.table("data_sync_logs").select("id", count="exact").limit(1).execute()
        print(f"SUCCESS: Table exists. Count: {res.count}")
    except Exception as e:
        print(f"FAILURE: Table check failed.")
        print(f"Error: {e}")
        if "relation" in str(e) and "does not exist" in str(e):
            print("CONCLUSION: Table 'data_sync_logs' is MISSING.")

if __name__ == "__main__":
    check_table()
