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

async def apply_migration():
    # Debug env
    print(f"Env keys: {[k for k in os.environ.keys() if 'URL' in k]}")
    
    db_url = getattr(settings, 'DATABASE_URL', None) or os.getenv('DATABASE_URL')
    if not db_url:
        print("Error: DATABASE_URL not found in settings.")
        # Try to construct from SUPABASE_URL if possible (unlikely without password)
        return

    sql_file = os.path.join(os.path.dirname(__file__), 'migrations', 'add_ai_refined_column.sql')
    with open(sql_file, 'r', encoding='utf-8') as f:
        sql = f.read()

    try:
        print(f"Connecting to database...")
        conn = await asyncpg.connect(db_url)
        print("Executing migration script...")
        await conn.execute(sql)
        print("✅ Migration completed successfully.")
        await conn.close()
    except Exception as e:
        print(f"❌ Migration failed: {e}")

if __name__ == "__main__":
    asyncio.run(apply_migration())
