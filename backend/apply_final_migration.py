import asyncio
import os
from app.core.supabase import supabase

async def run_final_migration():
    print("🚀 Running final schema migration via RPC...")
    
    # Test if exec_sql exists
    try:
        supabase.rpc('exec_sql', {'sql': 'SELECT 1'}).execute()
        print("✅ exec_sql RPC function is available.")
    except Exception as e:
        print(f"❌ exec_sql RPC function not found or failed: {e}")
        print("Please run the SQL manually in Supabase SQL Editor.")
        return

    sql_file = os.path.join(os.path.dirname(__file__), 'migrations', '20260126_final_schema_update.sql')
    if not os.path.exists(sql_file):
        print(f"❌ Migration file not found: {sql_file}")
        return

    with open(sql_file, 'r', encoding='utf-8') as f:
        sql = f.read()

    try:
        print(f"Executing migration: {os.path.basename(sql_file)}...")
        res = supabase.rpc('exec_sql', {'sql': sql}).execute()
        print("✅ Final migration via RPC successful!")
    except Exception as e:
        print(f"❌ Final migration via RPC failed: {e}")
        print("\n--- SQL CONTENT FOR MANUAL EXECUTION ---")
        print(sql)
        print("--- END SQL CONTENT ---")

if __name__ == "__main__":
    asyncio.run(run_final_migration())
