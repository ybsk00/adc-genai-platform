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
    # Try to get DATABASE_URL from env directly since it's not in Settings
    db_url = os.getenv("DATABASE_URL")
    if not db_url:
        # Fallback: check if settings has it (unlikely given previous error, but for safety)
        if hasattr(settings, "DATABASE_URL"):
            db_url = settings.DATABASE_URL
            
    if not db_url:
        print("Error: DATABASE_URL not found in environment variables or settings.")
        print("Available keys in env (partial):", [k for k in os.environ.keys() if "URL" in k or "DB" in k])
        return

    sql_file = os.path.join(os.path.dirname(__file__), 'migrations', '20260126_add_advanced_columns.sql')
    if not os.path.exists(sql_file):
        print(f"Error: Migration file not found at {sql_file}")
        return
        
    with open(sql_file, 'r', encoding='utf-8') as f:
        sql = f.read()

    try:
        print(f"Connecting to database...")
        conn = await asyncpg.connect(db_url)
        print(f"Executing migration script: {os.path.basename(sql_file)}")
        await conn.execute(sql)
        print("✅ Migration completed successfully.")
        await conn.close()
    except Exception as e:
        print(f"❌ Migration failed: {e}")

if __name__ == "__main__":
    asyncio.run(migrate())
