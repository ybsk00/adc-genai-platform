import os
import asyncio
import asyncpg
from dotenv import load_dotenv

# Load .env explicitly
env_path = os.path.join(os.path.dirname(__file__), '.env')
load_dotenv(env_path)

# Force reload config
from app.core import config
from importlib import reload
reload(config)
from app.core.config import settings

async def migrate():
    db_url = settings.DATABASE_URL
    if not db_url:
        # Try to construct from components if DATABASE_URL is missing
        # (Assuming SUPABASE_DB_URL or similar might be there, but let's check settings)
        print("Error: DATABASE_URL not found in settings.")
        return

    sql_file = os.path.join(os.path.dirname(__file__), 'db', 'migrate_to_gemini_768.sql')
    with open(sql_file, 'r', encoding='utf-8') as f:
        sql = f.read()

    try:
        print(f"Connecting to database...")
        conn = await asyncpg.connect(db_url)
        print("Executing migration script...")
        # asyncpg doesn't support multiple statements in one execute() easily if they are complex,
        # but for this script it should be fine.
        await conn.execute(sql)
        print("✅ Migration completed successfully.")
        await conn.close()
    except Exception as e:
        print(f"❌ Migration failed: {e}")

if __name__ == "__main__":
    asyncio.run(migrate())
