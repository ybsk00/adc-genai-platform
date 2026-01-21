import asyncio
import os
import sys
from dotenv import load_dotenv

# 1. Load .env manually
env_path = os.path.join(os.getcwd(), "backend", ".env")
if os.path.exists(env_path):
    with open(env_path, "r") as f:
        for line in f:
            if "=" in line and not line.startswith("#"):
                key, value = line.strip().split("=", 1)
                os.environ[key] = value.strip("'").strip('"')

# 2. Add backend to sys.path
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.services.ai_refiner import ai_refiner

async def verify_fallback():
    print("--- 🧪 AI Refiner Fallback Verification ---")
    
    # Test Drug: A complex ADC payload or a less common drug to trigger AI generation
    # Or we can just call generate_smiles_with_ai directly.
    test_drug = "Exatecan derivative (DXd)"
    print(f"\n[1] Testing AI SMILES Generation for: {test_drug}")
    
    try:
        result = await ai_refiner.generate_smiles_with_ai(test_drug)
        
        if result and "smiles_code" in result:
            print(f"   ✅ AI Generation Success!")
            print(f"      SMILES: {result['smiles_code']}")
            print(f"      Canonical: {result.get('canonical_smiles')}")
            print(f"      Source: {result.get('enrichment_source')}")
            
            # Double check with RDKit here just to be sure
            from rdkit import Chem
            mol = Chem.MolFromSmiles(result['smiles_code'])
            if mol:
                print("      ✅ RDKit Validation Passed (External Check)")
            else:
                print("      ❌ RDKit Validation Failed (External Check)")
                
        else:
            print(f"   ❌ AI Generation Failed: {result}")

    except Exception as e:
        print(f"   ❌ Test Error: {e}")

if __name__ == "__main__":
    asyncio.run(verify_fallback())
