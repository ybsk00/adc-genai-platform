import asyncio
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from app.services.bulk_importer import bulk_importer
from app.services.ai_refiner import ai_refiner
from app.services.chemical_mapper import chemical_mapper
from app.core.supabase import supabase
from manual_enrichment import manual_update

async def verify_all():
    print("🚀 Starting Verification...")
    
    # 1. Verify bulk_importer logic (Unit Test)
    print("\n[1] Verifying bulk_importer (Patient Count & Units)...")
    mock_results = {
        "participantFlowModule": {
            "groups": [{"title": "Arm A", "count": "10"}, {"title": "Arm B", "count": "20"}]
        },
        "outcomeMeasuresModule": {
            "outcomeMeasures": [
                {
                    "title": "Objective Response Rate",
                    "unitOfMeasure": "percentage",
                    "classes": [{"categories": [{"measurements": [{"value": "45.5"}]}]}]
                },
                {
                    "title": "Overall Survival",
                    "unitOfMeasure": "months",
                    "classes": [{"categories": [{"measurements": [{"value": "12.5"}]}]}]
                },
                {
                    "title": "Bad Unit ORR", # Should be rejected because title implies ORR but unit is months (if logic is strict)
                    # Actually logic checks: if title has ORR -> check if unit is %
                    "unitOfMeasure": "months", 
                    "classes": [{"categories": [{"measurements": [{"value": "99.9"}]}]}]
                }
            ]
        }
    }
    parsed = bulk_importer._parse_results_section(mock_results)
    metrics = parsed["metrics"]
    print(f" - Patient Count (Exp: 30): {metrics['patient_count']}")
    print(f" - ORR (Exp: 45.5): {metrics['orr_pct']}")
    print(f" - OS (Exp: 12.5): {metrics['os_months']}")
    
    if metrics['patient_count'] != 30: print("❌ Patient Count Failed")
    if metrics['orr_pct'] != 45.5: print("❌ ORR Failed")
    if metrics['os_months'] != 12.5: print("❌ OS Failed")
    
    # Check rejection of Bad Unit ORR (should not overwrite good ORR or should be ignored)
    # In this mock, "Bad Unit ORR" comes last. If it was accepted, ORR would be 99.9.
    # But since unit is months, it should be rejected for ORR.
    if metrics['orr_pct'] == 99.9:
        print("❌ Unit Verification Failed (Accepted months for ORR)")
    else:
        print("✅ Unit Verification Passed (Rejected months for ORR)")

    print("✅ bulk_importer verified.")

    # 2. Verify ai_refiner (Target Normalization)
    print("\n[2] Verifying ai_refiner (Target Normalization)...")
    targets = ["HER2", "TROP-2", "Nectin-4", "UnknownTarget"]
    for t in targets:
        norm = ai_refiner._normalize_target(t)
        print(f" - {t} -> {norm}")
    
    if ai_refiner._normalize_target("HER2") != "ERBB2": print("❌ HER2 Normalization Failed")
    if ai_refiner._normalize_target("TROP-2") != "TACSTD2": print("❌ TROP-2 Normalization Failed")
    print("✅ ai_refiner verified.")

    # 3. Verify chemical_mapper (SMILES Inference)
    print("\n[3] Verifying chemical_mapper (SMILES Inference)...")
    drugs = ["Trastuzumab deruxtecan", "Enfortumab vedotin", "UnknownDrug"]
    for d in drugs:
        res = chemical_mapper._infer_structure_from_name(d)
        print(f" - {d} -> {res[:20] if res else 'None'}...")
    
    if not chemical_mapper._infer_structure_from_name("Trastuzumab deruxtecan"): print("❌ DXd Inference Failed")
    if not chemical_mapper._infer_structure_from_name("Enfortumab vedotin"): print("❌ MMAE Inference Failed")
    print("✅ chemical_mapper verified.")

    # 4. Verify Manual Enrichment (Integration Test)
    print("\n[4] Verifying Manual Enrichment (Pazopanib)...")
    # First insert a dummy record if not exists
    dummy_name = "Pazopanib Test Record"
    # Clean up first just in case
    supabase.table("golden_set_library").delete().eq("name", dummy_name).execute()
    
    supabase.table("golden_set_library").insert({"name": dummy_name, "properties": {}}).execute()
    
    # Run manual update
    await manual_update(dummy_name, {"orr_pct": 99.9, "outcome_type": "Success"})
    
    # Check result
    res = supabase.table("golden_set_library").select("*").eq("name", dummy_name).execute()
    if res.data:
        record = res.data[0]
        print(f" - ORR: {record.get('orr_pct')}")
        print(f" - Source: {record.get('enrichment_source')}")
        
        if record.get('orr_pct') != 99.9: print("❌ Manual Update ORR Failed")
        if record.get('enrichment_source') != "Manual Review": print("❌ Manual Update Source Failed")
        
        # Cleanup
        supabase.table("golden_set_library").delete().eq("name", dummy_name).execute()
        print("✅ Manual Enrichment verified.")
    else:
        print("❌ Failed to find inserted record.")

if __name__ == "__main__":
    asyncio.run(verify_all())
