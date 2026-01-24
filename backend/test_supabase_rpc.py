import asyncio
from app.core.supabase import supabase

async def test_rpc():
    try:
        # Try to call a common helper function if it exists, or just a simple query
        # Since we can't easily run arbitrary SQL via the client without a helper function,
        # let's try to see if we can use the 'match_golden_set' which we know exists (or existed).
        
        # But our goal is to RUN the migration SQL.
        # If the user has enabled the 'exec_sql' RPC (common in some Supabase setups), we can use it.
        
        sql_file = 'backend/db/migrate_to_gemini_768.sql'
        with open(sql_file, 'r', encoding='utf-8') as f:
            sql = f.read()

        # Split SQL into individual statements as some clients/functions might not handle multiple
        statements = [s.strip() for s in sql.split(';') if s.strip()]
        
        print(f"Attempting to execute {len(statements)} statements via RPC...")
        
        # Note: This is a long shot if 'exec_sql' isn't defined.
        # Alternatively, we can try to use a python library that doesn't require psql but can connect via URL.
        # We already tried asyncpg and it failed because DATABASE_URL was missing from settings.
        
        # Let's try to find if there's any other way to get the connection string.
        pass

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    asyncio.run(test_rpc())
