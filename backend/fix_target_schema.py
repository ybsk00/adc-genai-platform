import os
import asyncio
import asyncpg
from dotenv import load_dotenv

# Load .env explicitly
env_path = os.path.join(os.path.dirname(__file__), '.env')
load_dotenv(env_path)

async def fix_schema():
    db_url = os.getenv('DATABASE_URL')
    if not db_url:
        # Fallback: Try to construct from SUPABASE credentials if DATABASE_URL is missing
        # This is a bit risky if connection pooling is involved, but worth a try for a fix script
        print("Warning: DATABASE_URL not found in .env. Checking for direct connection string...")
        return

    sql = """
    ALTER TABLE target_master 
    ADD COLUMN IF NOT EXISTS updated_at TIMESTAMPTZ DEFAULT NOW();

    ALTER TABLE antibody_library 
    ADD COLUMN IF NOT EXISTS updated_at TIMESTAMPTZ DEFAULT NOW();
    
    -- Add trigger for target_master
    DROP TRIGGER IF EXISTS handle_updated_at_target ON target_master;
    CREATE OR REPLACE FUNCTION update_updated_at_column()
    RETURNS TRIGGER AS $$
    BEGIN
        NEW.updated_at = NOW();
        RETURN NEW;
    END;
    $$ language 'plpgsql';
    
    CREATE TRIGGER handle_updated_at_target
    BEFORE UPDATE ON target_master
    FOR EACH ROW
    EXECUTE PROCEDURE update_updated_at_column();
    """

    try:
        print(f"Connecting to database...")
        conn = await asyncpg.connect(db_url)
        print("Executing schema fix...")
        await conn.execute(sql)
        print("✅ Schema fixed successfully: Added 'updated_at' to target_master and antibody_library.")
        await conn.close()
    except Exception as e:
        print(f"❌ Schema fix failed: {e}")

if __name__ == "__main__":
    asyncio.run(fix_schema())
