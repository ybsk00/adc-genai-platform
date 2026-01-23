import os
import logging
import asyncio
import asyncpg
from dotenv import load_dotenv

# Load .env explicitly
env_path = os.path.join(os.path.dirname(__file__), '.env')
print(f"Loading .env from: {env_path}")
load_dotenv(env_path)

# Force reload config to pick up env vars
from app.core import config
from importlib import reload
reload(config)
from app.core.config import settings

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def fix_embedding_dim():
    db_url = settings.DATABASE_URL
    if not db_url:
        logger.error("❌ DATABASE_URL not found in settings!")
        return

    try:
        logger.info("Connecting to database...")
        conn = await asyncpg.connect(db_url)
        
        # 1. Check if column exists and drop it (simplest way to change dim if data is not critical or incompatible)
        # The user said "DB가 입구를 닫아버린 상황" so likely no valid data with 1536 dims exists for this crawler yet, 
        # or we can just drop and re-add.
        logger.info("🗑️ Dropping existing embedding column (1536 dim)...")
        await conn.execute("ALTER TABLE commercial_reagents DROP COLUMN IF EXISTS embedding;")
        
        # 2. Add new column with 768 dims
        logger.info("✨ Adding new embedding column (768 dim)...")
        await conn.execute("ALTER TABLE commercial_reagents ADD COLUMN embedding vector(768);")
        
        # 3. Create index (optional but good for performance)
        logger.info("🔍 Creating hnsw index...")
        # Note: vector_cosine_ops is standard for cosine similarity
        await conn.execute("CREATE INDEX IF NOT EXISTS idx_commercial_reagents_embedding ON commercial_reagents USING hnsw (embedding vector_cosine_ops);")
        
        logger.info("✅ Successfully updated embedding dimension to 768.")
        await conn.close()
        
    except Exception as e:
        logger.error(f"❌ Migration failed: {e}")

if __name__ == "__main__":
    asyncio.run(fix_embedding_dim())
