import os
import asyncio
import asyncpg
import sys
from dotenv import load_dotenv

# Add backend directory to path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Load .env explicitly
env_path = os.path.join(os.path.dirname(__file__), '.env')
load_dotenv(env_path)

from app.core.config import settings

async def apply_migration():
    print(f"DEBUG: Env keys with URL: {[k for k in os.environ.keys() if 'URL' in k]}")
    db_url = getattr(settings, 'DATABASE_URL', None) or os.getenv('DATABASE_URL')
    if not db_url:
        print("❌ Error: DATABASE_URL not found in settings or .env")
        return

    sql_file = "add_enrichment_columns.sql"
    if not os.path.exists(sql_file):
        print(f"❌ Error: {sql_file} not found.")
        return

    with open(sql_file, 'r', encoding='utf-8') as f:
        sql = f.read()

    try:
        print(f"🔌 Connecting to database...")
        conn = await asyncpg.connect(db_url)
        print("🚀 Executing migration script...")
        await conn.execute(sql)
        print("✅ Migration completed successfully.")
        await conn.close()
    except Exception as e:
        print(f"❌ Migration failed: {e}")

if __name__ == "__main__":
    asyncio.run(apply_migration())
