import asyncio
import os
from app.core.supabase import supabase

async def test_exec_sql():
    print("Checking for 'exec_sql' RPC function...")
    try:
        # Simple test query
        res = supabase.rpc('exec_sql', {'sql': 'SELECT 1'}).execute()
        print(f"✅ exec_sql exists! Result: {res.data}")
        return True
    except Exception as e:
        print(f"❌ exec_sql not found or failed: {e}")
        return False

async def run_migration_via_rpc():
    if await test_exec_sql():
        print("Running migration via exec_sql...")
        sql_file = os.path.join(os.path.dirname(__file__), 'db', 'migrate_to_gemini_768.sql')
        with open(sql_file, 'r', encoding='utf-8') as f:
            sql = f.read()
        
        try:
            res = supabase.rpc('exec_sql', {'sql': sql}).execute()
            print("✅ Migration via RPC successful!")
        except Exception as e:
            print(f"❌ Migration via RPC failed: {e}")
    else:
        print("Cannot run migration via RPC. Please run the SQL manually in Supabase SQL Editor.")
        print("\n--- SQL CONTENT START ---")
        sql_file = os.path.join(os.path.dirname(__file__), 'db', 'migrate_to_gemini_768.sql')
        with open(sql_file, 'r', encoding='utf-8') as f:
            print(f.read())
        print("--- SQL CONTENT END ---")

if __name__ == "__main__":
    asyncio.run(run_migration_via_rpc())
