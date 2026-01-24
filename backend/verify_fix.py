import asyncio
import os
import sys
import logging
from dotenv import load_dotenv

# 환경 변수 로드 및 경로 설정
env_path = os.path.join(os.getcwd(), "backend", ".env")
load_dotenv(env_path)
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.services.rag_service import rag_service
from app.services.ai_refiner import ai_refiner
from app.core.supabase import supabase

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def verify_dimension():
    print("\n--- [1] Embedding Dimension Check ---")
    test_text = "Trastuzumab deruxtecan is an antibody-drug conjugate."
    embedding = await rag_service.generate_embedding(test_text)
    print(f"   Generated Embedding Length: {len(embedding)}")
    if len(embedding) == 1536:
        print("   ✅ Success: Dimension is 1536 (OpenAI Standard)")
    else:
        print(f"   ❌ Error: Dimension is {len(embedding)}, expected 1536")

async def verify_ai_refiner_commercial():
    print("\n--- [2] AI Refiner Commercial Reagent Support Check ---")
    mock_record = {
        "id": "test-id-123",
        "product_name": "HER2-directed ADC Toxin",
        "category": "ADC Toxin",
        "cas_number": "12345-67-8",
        "summary": "This is a high purity ADC toxin targeting HER2 antigen.",
        "source_name": "Ambeed",
        "properties": {"purity": "98%"}
    }
    
    print(f"   Testing Refiner with source_name: {mock_record['source_name']}")
    analysis = await ai_refiner.refine_single_record(mock_record)
    
    if analysis and "error" not in analysis:
        print("   ✅ Success: AI Refiner processed commercial reagent")
        print(f"      Extracted Target: {analysis.get('target')}")
        print(f"      Extracted Category: {analysis.get('category')}")
        print(f"      Relevance Score: {analysis.get('relevance_score')}")
    else:
        print(f"   ❌ Error: AI Refiner failed. Result: {analysis}")

async def main():
    await verify_dimension()
    await verify_ai_refiner_commercial()

if __name__ == "__main__":
    asyncio.run(main())
