import asyncio
import logging
import json
from app.core.supabase import supabase

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("FIX_AND_VERIFY")

async def fix_and_verify():
    print("\n" + "="*80)
    print("🛠️ [FIX] Updating Enrichment Source & Verifying Text Length")
    print("="*80 + "\n")

    target_drugs = ['Mylotarg', 'ELAHERE', 'ZYNLONTA', 'Besponsa', 'Blenrep']

    # 1. Force Update Source to 'open_fda_api' and Reset AI Status
    print(f"🔄 Updating {len(target_drugs)} records to source='open_fda_api'...")
    
    res = supabase.table("golden_set_library")\
        .update({
            "enrichment_source": "open_fda_api",
            "ai_refined": False,
            "relevance_score": 0,
            "outcome_type": "Unknown"
        })\
        .in_("name", target_drugs)\
        .execute()
    
    print(f"   ✅ Updated {len(res.data)} records.")

    # 2. Verify Text Lengths
    print("\n🔍 Verifying Text Retention in DB...")
    
    # Fetch updated records
    res = supabase.table("golden_set_library")\
        .select("name, properties, enrichment_source")\
        .in_("name", target_drugs)\
        .execute()
        
    # Write results to file to avoid encoding issues
    with open("fix_verify_results.txt", "w", encoding="utf-8") as f:
        f.write(f"\n{'Drug Name':<15} | {'Source':<15} | {'Indication Len':<15} | {'MoA Len':<15}\n")
        f.write("-" * 70 + "\n")
        
        for record in res.data:
            name = record['name']
            source = record['enrichment_source']
            props = record.get('properties', {})
            
            # Check both top-level and nested fda_label (just in case)
            ind = props.get('indication') or props.get('fda_label', {}).get('indication') or ""
            moa = props.get('mechanism_of_action') or props.get('fda_label', {}).get('mechanism_of_action') or ""
            
            # Handle list vs string (OpenFDA raw is list, our extract might be string)
            if isinstance(ind, list): ind = ind[0]
            if isinstance(moa, list): moa = moa[0]
            
            ind_len = len(ind)
            moa_len = len(moa)
            
            f.write(f"{name:<15} | {source:<15} | {ind_len:<15} | {moa_len:<15}\n")
            
            if ind_len < 1000:
                f.write(f"   ⚠️ WARNING: Indication text for {name} is short ({ind_len} chars).\n")

    print("Verification complete. Results written to fix_verify_results.txt")

    print("\n" + "="*80)
    print("✅ Ready for 'Run Full'. AI will now see 'open_fda_api' and full text.")

if __name__ == "__main__":
    asyncio.run(fix_and_verify())
