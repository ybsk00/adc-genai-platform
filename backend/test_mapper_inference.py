"""
Test Chemical Mapper Inference
- Product Name -> SMILES (Inference)
- Payload Alias -> SMILES
"""
import asyncio
import logging
from app.services.chemical_mapper import chemical_mapper

# Setup Logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def test_mapper():
    # 1. Test Inference from Name (Deruxtecan -> DXd)
    drug_name = "Trastuzumab deruxtecan"
    logger.info(f"🧪 Testing Inference for: {drug_name}")
    
    res = await chemical_mapper.enrich_with_commercial_data(drug_name)
    logger.info(f"Result: {res}")
    
    assert res.get("payload_smiles") is not None, "Payload SMILES inference failed"
    assert "Inference" in res.get("enrichment_source"), "Source should be Inference"

    # 2. Test Inference from Name (Vedotin -> MMAE)
    drug_name_2 = "Enfortumab vedotin"
    logger.info(f"\n🧪 Testing Inference for: {drug_name_2}")
    
    res_2 = await chemical_mapper.enrich_with_commercial_data(drug_name_2)
    logger.info(f"Result: {res_2}")
    
    assert res_2.get("payload_smiles") is not None, "Payload SMILES inference failed"
    
    logger.info("\n✅ All Mapper Tests Passed!")

if __name__ == "__main__":
    asyncio.run(test_mapper())
