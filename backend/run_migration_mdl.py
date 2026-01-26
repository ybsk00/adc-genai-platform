import os
import asyncio
import asyncpg
from dotenv import load_dotenv

# Load env
load_dotenv()

async def apply_migration():
    # Retrieve DATABASE_URL
    db_url = os.getenv('DATABASE_URL')
    if not db_url:
        # Fallback: Try to construct or warn
        print("❌ DATABASE_URL is missing! Cannot connect to DB directly.")
        # Sometimes connection string is different. Let's try to read it from other files if needed.
        return

    sql = """
    ALTER TABLE commercial_reagents ADD COLUMN IF NOT EXISTS mdl_number TEXT;
    """

    try:
        print(f"Connecting to DB...")
        conn = await asyncpg.connect(db_url)
        print("Executing: ADD COLUMN mdl_number")
        await conn.execute(sql)
        print("✅ Migration SUCCESS: mdl_number column added.")
        await conn.close()
    except Exception as e:
        print(f"❌ Migration FAILED: {e}")

if __name__ == "__main__":
    asyncio.run(apply_migration())
