import os
from app.core.supabase import supabase
import asyncio

async def apply_migration():
    print("Applying migration: add_enrichment_columns.sql")
    with open("add_enrichment_columns.sql", "r", encoding="utf-8") as f:
        sql = f.read()
    
    # Split by statement to execute one by one (Supabase RPC might handle multiple, but safer to split if using raw SQL via RPC)
    # However, supabase-py doesn't support raw SQL directly unless we have a 'exec_sql' RPC function.
    # Assuming we have one or we can use the 'rpc' method if a function exists.
    # If not, we might need to use a different approach.
    # Checking previous files, there seems to be 'apply_migration_fix.py' which uses a specific method.
    # Let's try to use the 'exec_sql' RPC if it exists, or just print instructions if we can't.
    
    # Actually, let's try to use the 'rpc' method to call 'exec_sql' which is a common pattern.
    try:
        # Try to execute as a single block if possible, or line by line
        # But wait, does the user have an 'exec_sql' function?
        # Let's assume yes or try to create it? No, I can't create it without SQL access.
        # I'll try to use the 'postgres' connection if I had it, but I only have supabase client.
        
        # Alternative: Use the 'apply_migration_advanced.py' pattern if it exists.
        # Let's just try to run it via a known RPC or just inform the user.
        # But I need to do it NOW.
        
        # Let's try to use the 'rpc' call 'execute_sql' or similar.
        res = supabase.rpc("exec_sql", {"sql_query": sql}).execute()
        print("Migration applied successfully via RPC.")
    except Exception as e:
        print(f"RPC failed: {e}")
        print("Trying alternative method (if available)...")
        # If RPC fails, we might not have permissions or the function.
        # In that case, I will ask the user to run it, but I should try my best first.

if __name__ == "__main__":
    asyncio.run(apply_migration())
