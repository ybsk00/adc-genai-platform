import os
import asyncio
import logging
from dotenv import load_dotenv
from app.core.supabase import supabase
from app.services.rag_service import rag_service

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

load_dotenv(os.path.join(os.path.dirname(__file__), '.env'))

async def debug_reindex():
    # Check counts
    res = supabase.table("golden_set_library").select("id", count="exact").limit(1).execute()
    count = res.count
    logger.info(f"Golden Set Library Count: {count}")
    
    res = supabase.table("golden_set_embeddings").select("id", count="exact").limit(1).execute()
    emb_count = res.count
    logger.info(f"Golden Set Embeddings Count: {emb_count}")
    
    if count > 0 and emb_count == 0:
        logger.info("Trying to index one item...")
        res = supabase.table("golden_set_library").select("*").limit(1).execute()
        item = res.data[0]
        logger.info(f"Item: {item['id']}")
        
        content = f"{item.get('description', '')}"
        embedding = await rag_service.generate_embedding(content)
        logger.info(f"Generated embedding length: {len(embedding)}")
        
        if embedding:
            data = {
                "source_id": item['id'],
                "chunk_content": content,
                "embedding": embedding
            }
            try:
                supabase.table("golden_set_embeddings").insert(data).execute()
                logger.info("✅ Inserted successfully.")
            except Exception as e:
                logger.error(f"❌ Insert failed: {e}")

if __name__ == "__main__":
    asyncio.run(debug_reindex())
