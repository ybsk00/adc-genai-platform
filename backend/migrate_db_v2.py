import os
import asyncio
import asyncpg
from dotenv import load_dotenv
import re

# Load .env explicitly
env_path = os.path.join(os.path.dirname(__file__), '.env')
load_dotenv(env_path)

# Force reload config
from app.core import config
from importlib import reload
reload(config)
from app.core.config import settings

async def migrate():
    # Try to construct DB URL from SUPABASE_URL if DATABASE_URL is missing
    # SUPABASE_URL: https://[project-ref].supabase.co
    # DB_URL: postgresql://postgres:[password]@db.[project-ref].supabase.co:5432/postgres
    
    url = settings.SUPABASE_URL
    project_ref = ""
    if url:
        match = re.search(r'https://(.*?)\.supabase\.co', url)
        if match:
            project_ref = match.group(1)
    
    # We need the password. Let's check if it's in .env under a different name
    # Common names: DB_PASSWORD, SUPABASE_DB_PASSWORD
    db_pass = os.getenv("DB_PASSWORD") or os.getenv("SUPABASE_DB_PASSWORD") or "your-password"
    
    db_url = f"postgresql://postgres:{db_pass}@db.{project_ref}.supabase.co:5432/postgres"
    
    print(f"Project Ref: {project_ref}")
    print(f"Attempting connection to: db.{project_ref}.supabase.co")

    sql_file = os.path.join(os.path.dirname(__file__), 'db', 'migrate_to_gemini_768.sql')
    with open(sql_file, 'r', encoding='utf-8') as f:
        sql = f.read()

    try:
        conn = await asyncpg.connect(db_url)
        print("Executing migration script...")
        await conn.execute(sql)
        print("✅ Migration completed successfully.")
        await conn.close()
    except Exception as e:
        print(f"❌ Migration failed: {e}")
        print("\nTrying to find DATABASE_URL in .env again with different method...")
        # If this fails, we might need to ask the user for the DB password or URL.

if __name__ == "__main__":
    asyncio.run(migrate())
