import asyncio
import logging
import json
import sys
from app.services.ai_refiner import ai_refiner
from app.services.openfda_service import openfda_service

# Configure logging to stdout
logging.basicConfig(level=logging.INFO, format='%(message)s', stream=sys.stdout)
logger = logging.getLogger("VERIFICATION")

async def verify_production_readiness():
    print("\n" + "="*80)
    print("🚀 [VERIFICATION] Production Readiness Check")
    print("="*80 + "\n")

    # 1. Text Insertion Check (OpenFDA Fetch)
    print("🔍 [TEST 1] Checking OpenFDA Data Fetch & Text Retention...")
    # We'll fetch Mylotarg specifically to see if we get the full text
    # Note: We are mocking the fetch result structure to test the EXTRACTOR logic, 
    # because calling the real OpenFDA API might be rate limited or slow, but we want to test the CODE logic.
    # However, to be 100% sure, let's try to fetch ONE real record if possible, or use a faithful mock.
    # Let's use the actual fetch_all_approved_adcs for a targeted search if possible.
    
    # Actually, let's just test the `extract_golden_info` method with a long dummy text to prove it doesn't truncate.
    long_text = "A" * 5000
    mock_fda_label = {
        "openfda": {
            "brand_name": ["MYLOTARG"],
            "generic_name": ["GEMTUZUMAB OZOGAMICIN"],
            "manufacturer_name": ["Pfizer"]
        },
        "indications_and_usage": [long_text],
        "mechanism_of_action": [long_text],
        "description": [long_text],
        "boxed_warning": [long_text]
    }
    
    extracted = openfda_service.extract_golden_info(mock_fda_label)
    props = extracted["properties"]
    
    print(f"   - Input Text Length: {len(long_text)}")
    print(f"   - Extracted Indication Length: {len(props['indication'])}")
    print(f"   - Extracted MoA Length: {len(props['mechanism_of_action'])}")
    
    if len(props['indication']) == len(long_text) and len(props['mechanism_of_action']) == len(long_text):
        print("   ✅ SUCCESS: Full text retained (No truncation detected).")
    else:
        print(f"   ❌ FAILED: Text truncated! Indication: {len(props['indication'])}, MoA: {len(props['mechanism_of_action'])}")

    # 2. Mylotarg Generic Name Fallback Test
    print("\n🔍 [TEST 2] Mylotarg SMILES & Target Extraction (Generic Name Fallback)...")
    
    # Create a record that mimics what would be in the DB
    record = {
        "id": "test-mylotarg",
        "name": "Mylotarg",
        "enrichment_source": "open_fda_api",
        "properties": {
            "generic_name": "GEMTUZUMAB OZOGAMICIN",
            "indication": "Mylotarg is indicated for the treatment of newly-diagnosed CD33-positive acute myeloid leukemia (AML)...",
            "mechanism_of_action": "Gemtuzumab ozogamicin is a CD33-directed antibody-drug conjugate (ADC)..."
        }
    }
    
    print(f"   - Target Drug: {record['name']}")
    print(f"   - Generic Name: {record['properties']['generic_name']}")
    
    # Run Refiner Logic
    # We manually call the steps to show the logs clearly
    
    # A. AI Analysis
    print("   👉 Running AI Analysis...")
    analysis = await ai_refiner.refine_single_record(record)
    
    print(f"   - AI Detected Target: {analysis.get('target')}")
    print(f"   - AI Relevance Score: {analysis.get('relevance_score')}")
    
    # B. SMILES Lookup (The Critical Step)
    print("   👉 Running SMILES Lookup...")
    drug_name = analysis.get("drug_name") or record['name']
    generic_name = analysis.get("generic_name") or record['properties']['generic_name']
    
    pubchem_data = await ai_refiner.enrich_with_pubchem(drug_name, generic_name)
    
    smiles = pubchem_data.get("smiles_code")
    source = pubchem_data.get("enrichment_source")
    
    print(f"   - SMILES Result: {smiles[:50]}..." if smiles else "   - SMILES Result: None")
    print(f"   - Source: {source}")
    
    if smiles and "GEMTUZUMAB" in generic_name:
         print("   ✅ SUCCESS: SMILES generated for Mylotarg using Generic Name!")
    else:
         print("   ❌ FAILED: SMILES not generated.")

    print("\n" + "="*80)

if __name__ == "__main__":
    asyncio.run(verify_production_readiness())
