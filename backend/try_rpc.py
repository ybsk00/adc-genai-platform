import os
import sys
from dotenv import load_dotenv

# Load .env explicitly
load_dotenv(os.path.join(os.path.dirname(__file__), '.env'))

# Add backend directory to path
sys.path.append(os.path.join(os.path.dirname(__file__), 'app'))

from app.core.supabase import supabase

def try_rpc():
    print("🧪 Testing exec_sql RPC...")
    sql = "ALTER TABLE commercial_reagents ADD COLUMN IF NOT EXISTS ai_refined BOOLEAN DEFAULT FALSE;"
    
    try:
        # Try 'exec_sql' or 'execute_sql'
        res = supabase.rpc("exec_sql", {"sql_query": sql}).execute()
        print(f"✅ RPC Success: {res}")
    except Exception as e:
        print(f"❌ RPC Failed: {e}")
        
        try:
            # Try alternative parameter name
            res = supabase.rpc("exec_sql", {"query": sql}).execute()
            print(f"✅ RPC (query param) Success: {res}")
        except Exception as e2:
             print(f"❌ RPC (query param) Failed: {e2}")

if __name__ == "__main__":
    try_rpc()
