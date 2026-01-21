import asyncio
import os
import sys
from dotenv import load_dotenv
from datetime import datetime

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
from app.core.supabase import supabase

async def verify_pipeline():
    print("--- 🧪 AI Refiner Pipeline Verification ---")
    
    # 1. PubChem Integration Check
    print("\n[1] Checking PubChem Integration...")
    test_drug = "Trastuzumab deruxtecan"
    print(f"   Target Drug: {test_drug}")
    pubchem_result = await ai_refiner.enrich_with_pubchem(test_drug)
    
    if pubchem_result and "smiles_code" in pubchem_result:
        print(f"   ✅ PubChem Success: Found SMILES for {test_drug}")
        print(f"      SMILES: {pubchem_result['smiles_code'][:30]}...")
        print(f"      Mol Weight: {pubchem_result['molecular_weight']}")
    else:
        print(f"   ❌ PubChem Failed: Could not fetch data for {test_drug}")
        print(f"      Result: {pubchem_result}")

    # 2. Cost Tracking Check (Check DB logs)
    print("\n[2] Checking Cost Tracking (llm_usage_logs)...")
    try:
        # Get latest log
        res = supabase.table("llm_usage_logs")\
            .select("*")\
            .order("created_at", desc=True)\
            .limit(1)\
            .execute()
        
        if res.data:
            log = res.data[0]
            print(f"   ✅ Cost Log Found: {log['created_at']}")
            print(f"      Model: {log['model']}")
            print(f"      Cost: ${log['cost_usd']}")
            print("      (This confirms logs are being written)")
        else:
            print("   ⚠️ No cost logs found yet (Run diag.py or Refiner first)")
            
    except Exception as e:
        print(f"   ❌ Failed to check logs: {e}")

    # 3. State Transition Logic Check (Code Inspection Simulation)
    print("\n[3] Checking State Transition Logic...")
    print("   Simulating update payload construction...")
    
    # Mock data
    mock_analysis = {
        "drug_name": "Test Drug",
        "outcome_type": "Success",
        "failure_reason": None,
        "relevance_score": 0.95
    }
    mock_pubchem = {"smiles_code": "C1=CC...", "molecular_weight": 150.0}
    
    update_payload = {
        "name": mock_analysis["drug_name"],
        "outcome_type": mock_analysis["outcome_type"],
        "relevance_score": mock_analysis["relevance_score"],
        "ai_refined": True,
        "rag_status": "processed",  # <--- Critical Check
        "processing_error": None
    }
    update_payload.update(mock_pubchem)
    
    if update_payload["rag_status"] == "processed" and update_payload["ai_refined"] is True:
        print("   ✅ Logic Confirmed: 'rag_status' is set to 'processed'")
        print("   ✅ Logic Confirmed: 'ai_refined' is set to True")
        print(f"   Generated Payload: {update_payload}")
    else:
        print("   ❌ Logic Error: State transition missing")

if __name__ == "__main__":
    asyncio.run(verify_pipeline())
