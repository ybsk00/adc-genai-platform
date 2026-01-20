import asyncio
import os
from dotenv import load_dotenv

# Load env vars
load_dotenv()

from supabase import create_client, Client

url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_KEY")

if not url or not key:
    print("Error: Missing SUPABASE_URL or SUPABASE_SERVICE_KEY")
    exit(1)

supabase: Client = create_client(url, key)

def test_insert():
    print(f"Testing insertion to 'golden_set_library'...")
    print(f"URL: {url}")
    
    data = {
        "name": "Complex Test Drug",
        "category": "clinical_trial",
        "description": "Test description with [Tag]",
        "status": "draft",
        "properties": {
            "nct_id": "NCT12345678",
            "title": "Complex Title",
            "conditions": ["Cancer", "Tumor"],
            "ai_analysis": {"extracted_name": "TestADC", "outcome_type": "Success"},
            "nested": {"deep": {"value": 123}}
        },
        "enrichment_source": "test_script",
        "raw_data": {"protocolSection": {"identificationModule": {"nctId": "NCT12345678"}}},
        "outcome_type": "Success",
        "failure_reason": None
    }
    
    try:
        # Try to insert
        response = supabase.table("golden_set_library").insert(data).execute()
        print("Success! Data inserted.")
        
        # Test Simple Select
        print("Testing Simple Select...")
        simple = supabase.table("golden_set_library").select("id").limit(1).execute()
        print(f"Simple Select Result: {simple.data}")

        # Test JSON Select with .contains()
        print("Testing JSON Select with .contains()...")
        nct_id = "NCT12345678"
        # Try using contains method on properties (which has nct_id at top level)
        existing = supabase.table("golden_set_library").select("id").contains("properties", {"nct_id": nct_id}).execute()
        print(f"JSON Select Result: {existing.data}")
        
        # Clean up
        print("Cleaning up...")
        record_id = response.data[0]['id']
        supabase.table("golden_set_library").delete().eq("id", record_id).execute()
        print("Test record deleted.")
        
    except Exception as e:
        print("\n!!! INSERTION FAILED !!!")
        print(f"Error Type: {type(e).__name__}")
        print(f"Error Details: {e}")
        
        # Check if it's a PostgREST error
        if hasattr(e, 'code'):
            print(f"Error Code: {e.code}")
        if hasattr(e, 'details'):
            print(f"Error Message: {e.details}")
            
        print("\nPossible Causes:")
        print("1. Table 'golden_set_library' does not exist. (Did you run the SQL script?)")
        print("2. Schema cache is stale. (Go to Supabase > API > Docs to reload, or restart project)")
        print("3. RLS Policy blocks insertion. (Check 'Service Role Full Access' policy)")

if __name__ == "__main__":
    test_insert()
