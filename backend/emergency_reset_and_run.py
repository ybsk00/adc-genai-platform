import asyncio
import logging
import json
import time
from app.core.supabase import supabase
from app.services.ai_refiner import ai_refiner

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("EMERGENCY_RUN")

async def emergency_reset_and_run():
    print("\n" + "="*80)
    print("🚨 [EMERGENCY] Hard Reset & Force Run for 5 Key Drugs")
    print("="*80 + "\n")

    target_drugs = ['Mylotarg', 'ELAHERE', 'ZYNLONTA', 'Besponsa', 'Blenrep']

    # 1. Hard Reset (Remove ai_analysis from properties)
    print("🧹 Step 1: Hard Resetting Records...")
    
    res = supabase.table("golden_set_library")\
        .select("*")\
        .in_("name", target_drugs)\
        .execute()
        
    for record in res.data:
        props = record.get("properties", {})
        if "ai_analysis" in props:
            del props["ai_analysis"] # Remove the old analysis
            
        # Ensure enrichment_source is open_fda_api
        # (User said force it)
        
        update_payload = {
            "ai_refined": False,
            "relevance_score": 0,
            "outcome_type": "Unknown",
            "enrichment_source": "open_fda_api", # FORCE SOURCE
            "properties": props,
            "updated_at": "now()"
        }
        
        supabase.table("golden_set_library")\
            .update(update_payload)\
            .eq("id", record["id"])\
            .execute()
            
    print(f"   ✅ Reset {len(res.data)} records. Ready for AI.")

    # 2. Force Run AI Refiner
    print("\n🤖 Step 2: Running AI Refiner (Force Mode)...")
    
    # We will manually call refine_single_record and update for each to show logs
    # This simulates what process_pending_records does but we want to see it happen now.
    
    for drug_name in target_drugs:
        print(f"\n   👉 Processing {drug_name}...")
        start_time = time.time()
        
        # Fetch fresh record
        res = supabase.table("golden_set_library").select("*").eq("name", drug_name).execute()
        if not res.data:
            print(f"      ❌ Record not found!")
            continue
            
        record = res.data[0]
        
        # Verify text length before run
        props = record.get("properties", {})
        ind = props.get("indication") or props.get("fda_label", {}).get("indication") or ""
        if isinstance(ind, list): ind = ind[0]
        print(f"      - Input Text Length: {len(ind)} chars")
        
        # Run Analysis
        # Note: refine_single_record does NOT update DB. It returns the analysis.
        # We must update DB manually here to complete the cycle.
        
        analysis = await ai_refiner.refine_single_record(record)
        
        duration = time.time() - start_time
        print(f"      - Analysis Time: {duration:.2f}s")
        
        if analysis and "error" not in analysis:
            target = analysis.get("target")
            score = analysis.get("relevance_score")
            print(f"      - Extracted Target: {target}")
            print(f"      - Relevance Score: {score}")
            
            # Run SMILES enrichment (Manual call to match process_pending_records logic)
            # We use the updated logic we added to ai_refiner
            drug_name_extracted = analysis.get("drug_name") or record["name"]
            generic_name = analysis.get("generic_name") or props.get("generic_name")
            
            pubchem_data = await ai_refiner.enrich_with_pubchem(drug_name_extracted, generic_name)
            smiles = pubchem_data.get("smiles_code") if pubchem_data else None
            print(f"      - SMILES: {'✅ Found' if smiles else '❌ Failed'}")
            
            # Update DB
            updated_props = record.get("properties", {})
            updated_props["ai_analysis"] = analysis
            
            update_payload = {
                "ai_refined": True,
                "relevance_score": score,
                "outcome_type": analysis.get("outcome_type", "Unknown"),
                "properties": updated_props
            }
            if pubchem_data:
                update_payload.update(pubchem_data)
                
            supabase.table("golden_set_library").update(update_payload).eq("id", record["id"]).execute()
            print("      💾 DB Updated.")
            
        else:
            print(f"      ❌ Analysis Failed: {analysis}")

    print("\n" + "="*80)
    print("✅ Emergency Run Complete.")

if __name__ == "__main__":
    asyncio.run(emergency_reset_and_run())
