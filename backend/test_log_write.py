import os
from dotenv import load_dotenv
from supabase import create_client, Client

load_dotenv()

url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_KEY")

def test_log_write():
    print("=== Testing Write to 'data_sync_logs' ===")
    supabase: Client = create_client(url, key)
    
    try:
        # Try to insert a test log
        data = {
            "source_id": "test_script",
            "status": "running",
            "records_synced": 0,
            "records_drafted": 0,
            "error_message": "Test Log from Local Script"
        }
        res = supabase.table("data_sync_logs").insert(data).execute()
        
        if res.data:
            print("SUCCESS: Successfully wrote to 'data_sync_logs'.")
            print(f"Log ID: {res.data[0]['id']}")
            
            # Clean up
            supabase.table("data_sync_logs").delete().eq("id", res.data[0]['id']).execute()
            print("Test log deleted.")
        else:
            print("FAILURE: Insert returned no data.")
            
    except Exception as e:
        print(f"FAILURE: Write failed.")
        print(f"Error: {e}")
        if "permission denied" in str(e):
            print("CONCLUSION: RLS is blocking writes to 'data_sync_logs'.")

if __name__ == "__main__":
    test_log_write()
