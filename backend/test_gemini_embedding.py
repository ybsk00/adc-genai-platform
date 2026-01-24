import asyncio
import os
from dotenv import load_dotenv

# Load .env
env_path = os.path.join(os.path.dirname(__file__), '.env')
load_dotenv(env_path)

# Import rag_service
from app.services.rag_service import rag_service

async def test_embedding():
    print("Testing Gemini Embedding Generation...")
    text = "This is a test for ADC platform RAG dimension unification."
    embedding = await rag_service.generate_embedding(text)
    
    if embedding:
        print(f"✅ Success! Embedding generated.")
        print(f"📏 Dimension: {len(embedding)}")
        if len(embedding) == 768:
            print("✅ Dimension is correctly 768.")
        else:
            print(f"❌ Dimension mismatch: Expected 768, got {len(embedding)}")
    else:
        print("❌ Failed to generate embedding.")

if __name__ == "__main__":
    asyncio.run(test_embedding())
