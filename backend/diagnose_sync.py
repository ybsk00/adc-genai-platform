import os
import json
import base64
from dotenv import load_dotenv
from supabase import create_client, Client

# Load env vars
load_dotenv()

url = os.getenv("SUPABASE_URL")
key = os.getenv("SUPABASE_SERVICE_KEY")

output = []

def log(msg):
    print(msg)
    output.append(str(msg))

def decode_jwt_role(token):
    try:
        # JWT is header.payload.signature
        payload_part = token.split('.')[1]
        # Add padding if needed
        payload_part += '=' * (-len(payload_part) % 4)
        payload = json.loads(base64.b64decode(payload_part).decode('utf-8'))
        return payload.get('role', 'unknown')
    except Exception as e:
        return f"error_decoding: {e}"

def run_diagnostics():
    log("=== START DIAGNOSTICS ===")
    
    if not url or not key:
        log("FATAL: Missing SUPABASE_URL or SUPABASE_SERVICE_KEY")
        return

    # Check Key Role
    role = decode_jwt_role(key)
    log(f"API Key Role: {role}")
    if role != 'service_role':
        log("WARNING: You are NOT using the service_role key. INSERTs might fail if RLS only allows service_role.")

    try:
        supabase: Client = create_client(url, key)
        log("Supabase client initialized.")
    except Exception as e:
        log(f"FATAL: Client init failed: {e}")
        return

    # 1. Check Table Existence & Read Permission
    log("\n[Test 1] Checking 'golden_set_library' existence (SELECT)...")
    try:
        res = supabase.table("golden_set_library").select("id", count="exact").limit(1).execute()
        log(f"SUCCESS: Table exists. Count: {res.count}")
    except Exception as e:
        log(f"FAILURE: Select failed. Error: {e}")
        return

    # 2. Check Write Permission (Simple)
    log("\n[Test 2] Checking Simple INSERT...")
    simple_data = {
        "name": "Diag Simple",
        "description": "Simple test",
        "status": "draft"
    }
    simple_id = None
    try:
        res = supabase.table("golden_set_library").insert(simple_data).execute()
        if res.data:
            simple_id = res.data[0]['id']
            log(f"SUCCESS: Simple insert worked. ID: {simple_id}")
        else:
            log("FAILURE: Insert returned no data.")
    except Exception as e:
        log(f"FAILURE: Simple insert failed. Error: {e}")

    # 3. Check JSONB Write (Complex)
    log("\n[Test 3] Checking Complex JSONB INSERT...")
    complex_data = {
        "name": "Diag Complex",
        "description": "Complex test",
        "status": "draft",
        "properties": {"test_key": "test_value", "nested": {"a": 1}},
        "raw_data": {"source": "diagnostic"}
    }
    complex_id = None
    try:
        res = supabase.table("golden_set_library").insert(complex_data).execute()
        if res.data:
            complex_id = res.data[0]['id']
            log(f"SUCCESS: Complex insert worked. ID: {complex_id}")
        else:
            log("FAILURE: Complex insert returned no data.")
    except Exception as e:
        log(f"FAILURE: Complex insert failed. Error: {e}")

    # 4. Check Duplicate Logic (JSON Filter)
    log("\n[Test 4] Checking JSON Filter (Duplicate Check)...")
    try:
        # Test .contains
        res = supabase.table("golden_set_library").select("id").contains("properties", {"test_key": "test_value"}).execute()
        log(f"SUCCESS: .contains() query worked. Found: {len(res.data)}")
    except Exception as e:
        log(f"FAILURE: .contains() query failed. Error: {e}")

    # Cleanup
    log("\n[Cleanup] Deleting test records...")
    try:
        ids_to_delete = []
        if simple_id: ids_to_delete.append(simple_id)
        if complex_id: ids_to_delete.append(complex_id)
        
        if ids_to_delete:
            supabase.table("golden_set_library").delete().in_("id", ids_to_delete).execute()
            log("Cleanup successful.")
    except Exception as e:
        log(f"Cleanup warning: {e}")

    log("\n=== END DIAGNOSTICS ===")

if __name__ == "__main__":
    run_diagnostics()
    with open("backend/final_diag.txt", "w", encoding="utf-8") as f:
        f.write("\n".join(output))
