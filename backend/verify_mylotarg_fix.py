import asyncio
import logging
import json
from app.core.supabase import supabase
from app.services.ai_refiner import ai_refiner

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def verify_mylotarg():
    print("🔍 Verifying Mylotarg Fix...")
    
    # 1. Fetch Mylotarg record (or insert if missing for test)
    res = supabase.table("golden_set_library")\
        .select("*")\
        .ilike("name", "Mylotarg")\
        .execute()
    
    record = None
    if res.data:
        record = res.data[0]
        print(f"✅ Found existing record: {record['name']} (ID: {record['id']})")
    else:
        print("⚠️ Mylotarg not found in DB. Creating mock record...")
        # Mock record simulating OpenFDA data
        mock_props = {
            "fda_label": {
                "brand_name": ["MYLOTARG"],
                "generic_name": ["GEMTUZUMAB OZOGAMICIN"],
                "mechanism_of_action": ["MYLOTARG is a CD33-directed antibody-drug conjugate (ADC). The antibody portion (gemtuzumab) binds specifically to the CD33 antigen..."],
                "indications_and_usage": ["MYLOTARG is indicated for the treatment of newly-diagnosed CD33-positive acute myeloid leukemia (AML)..."],
                "description": ["MYLOTARG (gemtuzumab ozogamicin) is an antibody-drug conjugate..."]
            },
            "generic_name": "GEMTUZUMAB OZOGAMICIN"
        }
        record = {
            "name": "Mylotarg",
            "enrichment_source": "open_fda_api",
            "properties": mock_props,
            "id": "mock-id-123"
        }

    # 2. Run Refiner
    print("\n🤖 Running AI Refiner on record...")
    # We need to simulate the process_pending_records logic partially or just call refine_single_record
    # But refine_single_record doesn't do the DB update or SMILES lookup logic fully (that's in process_pending_records).
    # However, ai_refiner.process_pending_records is hard to run for just one specific item without modifying it.
    # So I will manually invoke the steps: refine -> enrich
    
    # Step 1: AI Analysis
    analysis = await ai_refiner.refine_single_record(record)
    print(f"\n📊 AI Analysis Result:\n{json.dumps(analysis, indent=2)}")
    
    drug_name = analysis.get("drug_name")
    generic_name = analysis.get("generic_name") or record.get("properties", {}).get("generic_name")
    
    print(f"\n🧪 Drug Name: {drug_name}")
    print(f"🧬 Generic Name: {generic_name}")
    
    # Step 2: SMILES Lookup (The core fix)
    print("\n🔬 Testing SMILES Lookup (with Fallback)...")
    pubchem_data = await ai_refiner.enrich_with_pubchem(drug_name, generic_name)
    
    print(f"\n📦 PubChem Result:\n{json.dumps(pubchem_data, indent=2)}")
    
    if pubchem_data and pubchem_data.get("smiles_code"):
        print("\n✅ SUCCESS: SMILES Generated!")
        if pubchem_data.get("enrichment_source") == "PubChem":
            print("   (Found via PubChem - likely using Generic Name fallback)")
    else:
        print("\n❌ FAILED: SMILES not found.")

if __name__ == "__main__":
    asyncio.run(verify_mylotarg())
