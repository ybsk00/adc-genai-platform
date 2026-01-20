import os
from dotenv import load_dotenv
from supabase import create_client, Client

load_dotenv()

url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_KEY")

def check_logs():
    print("=== Checking Data Sync Logs ===")
    supabase: Client = create_client(url, key)
    
    try:
        # Fetch latest 5 failed logs
        res = supabase.table("data_sync_logs") \
            .select("*") \
            .eq("status", "failed") \
            .order("started_at", desc=True) \
            .limit(5) \
            .execute()
            
        if not res.data:
            print("No failed logs found in 'data_sync_logs'.")
            print("Possibilities:")
            print("1. The backend code with DB logging hasn't been deployed/restarted.")
            print("2. The sync job hasn't run since the code update.")
        else:
            print(f"Found {len(res.data)} failed logs:")
            for log in res.data:
                print(f"--- Log ID: {log.get('id')} ---")
                print(f"Time: {log.get('started_at')}")
                print(f"Error: {log.get('error_message')}")
                print("-----------------------------")
                
    except Exception as e:
        print(f"Error fetching logs: {e}")

if __name__ == "__main__":
    check_logs()
