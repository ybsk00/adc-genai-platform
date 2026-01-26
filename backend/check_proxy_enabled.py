import os
from dotenv import load_dotenv

def check_proxy():
    env_path = os.path.join("backend", ".env")
    load_dotenv(env_path)
    proxy_enabled = os.getenv("PROXY_ENABLED")
    print(f"PROXY_ENABLED={proxy_enabled}")

if __name__ == "__main__":
    check_proxy()
