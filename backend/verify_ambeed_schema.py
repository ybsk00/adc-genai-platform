import asyncio
import logging
import sys
import os
from unittest.mock import MagicMock, patch

# Add backend directory to sys.path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Mock supabase before importing crawler to avoid real DB calls during schema check
with patch('app.services.ambeed_crawler.supabase') as mock_supabase:
    from app.services.ambeed_crawler import ambeed_crawler

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def verify_schema_match():
    print("🧪 Starting Ambeed Schema Verification...")

    # Mock data that would come from the page
    mock_raw_data = {
        "ambeed_cat_no": "TEST-123",
        "product_name": "Test ADC Toxin",
        "product_url": "http://test.com",
        "category": "ADC Toxins",
        "smiles_code": "C1=CC=C(C=C1)C(=O)O", # Simple SMILES
        "cas_number": "123-45-6",
        "body_text": "This is a test toxin. Payload is MMAE. Linker is Val-Cit.",
        "crawled_at": "2026-01-25T00:00:00",
        "properties": {}
    }

    # Mock Gemini response
    mock_ai_response = {
        "payload_smiles": "PAYLOAD_SMILES_TEST",
        "linker_smiles": "LINKER_SMILES_TEST",
        "full_smiles": "FULL_SMILES_TEST",
        "target": "Test Target",
        "summary": "Test Summary"
    }

    # Patch _enrich_with_gemini to return mock data
    with patch.object(ambeed_crawler, '_enrich_with_gemini', return_value=mock_ai_response):
        print("🤖 Mocking AI response...")
        
        # Run _enrich_and_prepare_item
        final_data = await ambeed_crawler._enrich_and_prepare_item(mock_raw_data)
        
        if not final_data:
            print("❌ _enrich_and_prepare_item returned None!")
            return

        print("\n📊 Generated Data Structure:")
        for k, v in final_data.items():
            print(f"  - {k}: {v}")

        # Verify Schema Keys
        required_keys = [
            "ambeed_cat_no", "product_name", "category", "cas_number", 
            "smiles_code", "payload_smiles", "linker_smiles", "full_smiles",
            "ai_refined", "target", "summary"
        ]
        
        missing_keys = [key for key in required_keys if key not in final_data]
        
        if missing_keys:
            print(f"\n❌ Missing Keys for Schema: {missing_keys}")
        else:
            print("\n✅ All required schema keys are present.")

        # Verify Values
        if final_data["payload_smiles"] == "PAYLOAD_SMILES_TEST" and \
           final_data["linker_smiles"] == "LINKER_SMILES_TEST" and \
           final_data["full_smiles"] == "FULL_SMILES_TEST":
            print("✅ 3-Part SMILES correctly mapped.")
        else:
            print("❌ 3-Part SMILES mapping failed.")

if __name__ == "__main__":
    asyncio.run(verify_schema_match())
