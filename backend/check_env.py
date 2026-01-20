import os

def check_env():
    env_path = ".env"
    if not os.path.exists(env_path):
        print(f"Error: {env_path} does not exist.")
        return

    print(f"Found {env_path}.")
    
    required_keys = ["SUPABASE_URL", "SUPABASE_SERVICE_KEY", "OPENAI_API_KEY"]
    found_keys = []
    
    with open(env_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                parts = line.split("=", 1)
                key = parts[0].strip()
                val = parts[1].strip()
                if key in required_keys and val:
                    found_keys.append(key)
    
    print("Keys found in .env (Non-empty):")
    for key in required_keys:
        if key in found_keys:
            print(f" - {key}: [PRESENT]")
        else:
            print(f" - {key}: [MISSING]")

if __name__ == "__main__":
    check_env()
