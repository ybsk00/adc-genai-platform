import os
import asyncio
import asyncpg
from dotenv import load_dotenv

# Load .env
env_path = os.path.join(os.path.dirname(__file__), '.env')
load_dotenv(env_path)

async def add_columns():
    db_url = os.getenv('DATABASE_URL')
    if not db_url:
        print("Error: DATABASE_URL not found in .env")
        return

    sql = """
    ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS binding_affinity TEXT;
    ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS isotype TEXT;
    ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS host_species TEXT;
    ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS orr_pct TEXT;
    ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS os_months TEXT;
    ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS pfs_months TEXT;
    ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS is_manual_override BOOLEAN DEFAULT FALSE;
    """

    try:
        print(f"Connecting to database...")
        conn = await asyncpg.connect(db_url)
        print("Executing column addition script...")
        await conn.execute(sql)
        print("✅ golden_set_library table updated successfully.")
        await conn.close()
    except Exception as e:
        print(f"❌ Column addition failed: {e}")

if __name__ == "__main__":
    asyncio.run(add_columns())
