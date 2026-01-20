import os
import httpx
import asyncio

async def test_curl():
    env_path = ".env"
    url = ""
    key = ""
    
    # Read .env manually
    if os.path.exists(env_path):
        with open(env_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip().startswith("SUPABASE_URL="):
                    url = line.split("=", 1)[1].strip()
                elif line.strip().startswith("SUPABASE_SERVICE_KEY="):
                    key = line.split("=", 1)[1].strip()
    
    if not url or not key:
        print("Error: Missing URL or Key in .env")
        return

    print(f"Testing URL: {url}")
    
    # Clean URL
    if url.endswith("/"):
        url = url[:-1]
        
    target_endpoint = f"{url}/rest/v1/golden_set_library"
    
    headers = {
        "apikey": key,
        "Authorization": f"Bearer {key}",
        "Content-Type": "application/json"
    }
    
    print(f"Requesting: GET {target_endpoint}")
    
    async with httpx.AsyncClient() as client:
        try:
            # Test GET
            response = await client.get(target_endpoint, headers=headers, params={"select": "*", "limit": 1})
            print(f"GET Status Code: {response.status_code}")
            print(f"GET Body: {response.text}")
            
            # Test POST (Insert)
            print(f"Requesting: POST {target_endpoint}")
            payload = {
                "name": "Curl Test Drug",
                "status": "draft"
            }
            response = await client.post(target_endpoint, headers=headers, json=payload)
            print(f"POST Status Code: {response.status_code}")
            print(f"POST Body: {response.text}")
            
        except Exception as e:
            print(f"Request Failed: {e}")

if __name__ == "__main__":
    asyncio.run(test_curl())
