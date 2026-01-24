import os
import logging
from dotenv import load_dotenv
from app.core.supabase import supabase

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

load_dotenv(os.path.join(os.path.dirname(__file__), '.env'))

def verify_dims():
    try:
        # Check commercial_reagents
        response = supabase.table("commercial_reagents").select("id, embedding").limit(5).execute()
        if response.data:
            for item in response.data:
                emb = item.get("embedding")
                if emb:
                    logger.info(f"✅ Commercial Item {item['id']} has embedding length: {len(emb)}")
        else:
            logger.info("⚠️ No data in commercial_reagents.")

        # Check golden_set_embeddings
        response = supabase.table("golden_set_embeddings").select("id, embedding").limit(5).execute()
        if response.data:
             for item in response.data:
                emb = item.get("embedding")
                if emb:
                    logger.info(f"✅ Golden Set Item {item['id']} has embedding length: {len(emb)}")
                    if len(emb) != 768:
                        logger.error(f"❌ Dimension mismatch! Expected 768, got {len(emb)}")
        else:
            logger.info("⚠️ No data in golden_set_embeddings.")
                
    except Exception as e:
        logger.error(f"❌ Verification failed: {e}")

if __name__ == "__main__":
    verify_dims()
