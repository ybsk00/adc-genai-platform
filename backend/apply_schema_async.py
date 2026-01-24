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

async def apply_schema():
    db_url = getattr(settings, 'DATABASE_URL', None) or os.getenv('DATABASE_URL')
    if not db_url:
        print("Error: DATABASE_URL not found in settings.")
        return

    sql_file = os.path.join(os.path.dirname(__file__), 'update_commercial_reagents_schema.sql')
    with open(sql_file, 'r', encoding='utf-8') as f:
        sql = f.read()

    try:
        print(f"Connecting to database...")
        conn = await asyncpg.connect(db_url)
        print("Executing schema update script...")
        await conn.execute(sql)
        print("✅ Schema update completed successfully.")
        await conn.close()
    except Exception as e:
        print(f"❌ Schema update failed: {e}")

if __name__ == "__main__":
    asyncio.run(apply_schema())
