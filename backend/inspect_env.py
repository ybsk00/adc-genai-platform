import os
from dotenv import load_dotenv

env_path = os.path.join(os.path.dirname(__file__), '.env')
print(f"Loading .env from: {env_path}")
load_dotenv(env_path)

print("Available keys in .env:")
for key in os.environ:
    if "URL" in key or "DB" in key or "KEY" in key:
        print(f" - {key}")
