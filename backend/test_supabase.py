from app.core.supabase import supabase
import uuid

def test_supabase_connection():
    print("--- Testing Supabase Connection ---")
    
    # 1. Test Select
    try:
        print("Attempting to select from golden_set_library...")
        response = supabase.table("golden_set_library").select("*").limit(1).execute()
        print(f"Select successful. Data: {response.data}")
    except Exception as e:
        print(f"Select Failed: {e}")
        return

    # 2. Test Insert (Draft)
    try:
        print("Attempting to insert a test record...")
        test_id = str(uuid.uuid4())
        test_record = {
            "name": f"Test Drug {test_id[:8]}",
            "category": "test",
            "description": "Connection test record",
            "status": "draft",
            "enrichment_source": "manual_test",
            "raw_data": {"test": True}
        }
        
        response = supabase.table("golden_set_library").insert(test_record).execute()
        print("Insert successful.")
        
        # Cleanup (Optional, or just leave it as draft)
        # supabase.table("golden_set_library").delete().eq("name", test_record["name"]).execute()
        # print("Cleanup successful.")
        
    except Exception as e:
        print(f"Insert Failed: {e}")
        if hasattr(e, 'code'):
            print(f"Error Code: {e.code}")
        if hasattr(e, 'details'):
            print(f"Error Details: {e.details}")
        if hasattr(e, 'message'):
            print(f"Error Message: {e.message}")

if __name__ == "__main__":
    from app.core.config import settings
    print(f"Target URL: {settings.SUPABASE_URL}")
    test_supabase_connection()
