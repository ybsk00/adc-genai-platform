import asyncio
import os
from dotenv import load_dotenv
from app.core.supabase import supabase
from app.services.rag_service import rag_service

# Load .env
env_path = os.path.join(os.path.dirname(__file__), '.env')
load_dotenv(env_path)

async def test_db_storage():
    print("Testing 768-dim Vector DB Storage...")
    text = "Test storage for 768-dim vector."
    embedding = await rag_service.generate_embedding(text)
    
    if not embedding or len(embedding) != 768:
        print(f"❌ Failed to generate 768-dim embedding. Got {len(embedding) if embedding else 0}")
        return

    try:
        # Try to insert into commercial_reagents (using a dummy record)
        test_data = {
            "ambeed_cat_no": "TEST-768-DIM",
            "product_name": "768 Dimension Test Product",
            "embedding": embedding,
            "source_name": "Test"
        }
        
        print("Attempting to upsert test record into commercial_reagents...")
        res = supabase.table("commercial_reagents").upsert(test_data, on_conflict="ambeed_cat_no").execute()
        
        if res.data:
            print("✅ Success! 768-dim vector stored in commercial_reagents.")
            # Clean up
            supabase.table("commercial_reagents").delete().eq("ambeed_cat_no", "TEST-768-DIM").execute()
            print("🧹 Test record cleaned up.")
        else:
            print(f"❌ Failed to store vector. Response: {res}")
            
    except Exception as e:
        print(f"❌ DB Storage Test failed: {e}")

if __name__ == "__main__":
    asyncio.run(test_db_storage())
