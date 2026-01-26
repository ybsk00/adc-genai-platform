import os
from dotenv import load_dotenv

def list_env_keys():
    env_path = os.path.join("backend", ".env")
    if not os.path.exists(env_path):
        print(f"{env_path} not found")
        return
    
    load_dotenv(env_path)
    print("Environment variable keys in .env:")
    with open(env_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key = line.split("=")[0]
                print(key)

if __name__ == "__main__":
    list_env_keys()
