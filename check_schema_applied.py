import os
import sys
from dotenv import load_dotenv

# 1. Load .env manually
env_path = os.path.join(os.getcwd(), "backend", ".env")
if os.path.exists(env_path):
    with open(env_path, "r") as f:
        for line in f:
            if "=" in line and not line.startswith("#"):
                key, value = line.strip().split("=", 1)
                os.environ[key] = value.strip("'").strip('"')
else:
    print("Error: .env not found")
    sys.exit(1)

# 2. Add backend to sys.path
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.core.supabase import supabase

def check_schema():
    print("--- Checking Schema ---")
    try:
        # Try to select the new columns from a single record
        res = supabase.table("golden_set_library")\
            .select("id, relevance_score, rag_status, molecular_weight, canonical_smiles")\
            .limit(1)\
            .execute()
        
        print("Schema check PASSED. Columns exist.")
        print(f"Sample data: {res.data}")
        
    except Exception as e:
        print(f"Schema check FAILED: {e}")
        print("Please run the SQL provided in the previous step in Supabase SQL Editor.")

if __name__ == "__main__":
    check_schema()
