import os
from dotenv import load_dotenv

load_dotenv()

print("Available Environment Variables:")
for key, value in os.environ.items():
    if "KEY" in key or "SECRET" in key or "PASSWORD" in key or "URL" in key or "DB" in key:
        print(f"{key}: {'*' * 8}") # Hide values
    else:
        print(key)

print("-" * 20)
print(f"DATABASE_URL exists: {'DATABASE_URL' in os.environ}")
if 'DATABASE_URL' in os.environ:
    print(f"DATABASE_URL length: {len(os.environ['DATABASE_URL'])}")
