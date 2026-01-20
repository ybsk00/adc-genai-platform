import os
import json
import base64

def check_key():
    env_path = ".env"
    if not os.path.exists(env_path):
        print("Error: .env not found")
        return

    key = ""
    with open(env_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip().startswith("SUPABASE_SERVICE_KEY="):
                key = line.split("=", 1)[1].strip()
                break
    
    if not key:
        print("Error: SUPABASE_SERVICE_KEY not found in .env")
        return

    print(f"Key found. Length: {len(key)}")
    
    try:
        # JWT is header.payload.signature
        parts = key.split(".")
        if len(parts) != 3:
            print("Error: Key does not look like a valid JWT (does not have 3 parts).")
            return
            
        # Decode payload (2nd part)
        payload_b64 = parts[1]
        # Add padding if needed
        payload_b64 += "=" * ((4 - len(payload_b64) % 4) % 4)
        
        payload_json = base64.urlsafe_b64decode(payload_b64).decode('utf-8')
        payload = json.loads(payload_json)
        
        print("--- Key Payload Analysis ---")
        print(f"Role: {payload.get('role')}")
        print(f"Iss: {payload.get('iss')}")
        
        if payload.get('role') == 'service_role':
            print("SUCCESS: This looks like a valid SERVICE_ROLE key.")
        elif payload.get('role') == 'anon':
            print("WARNING: This is an ANON (Public) key. You need the SERVICE_ROLE key.")
        else:
            print(f"WARNING: Unknown role '{payload.get('role')}'.")
            
    except Exception as e:
        print(f"Error decoding key: {e}")

if __name__ == "__main__":
    check_key()
