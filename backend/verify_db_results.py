import asyncio
import logging
import json
from app.core.supabase import supabase

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def verify_db_results():
    print("\n" + "="*80)
    print("📊 [VERIFICATION] Database Result Check")
    print("="*80 + "\n")

    # List of drugs to check (from user's previous failure list)
    target_drugs = [
        "Mylotarg", 
        "HERZUMA", 
        "Besponsa", 
        "ELAHERE", 
        "ZYNLONTA", 
        "Ontruzant", 
        "Blenrep"
    ]
    
    # Fetch records
    print(f"🔍 Checking {len(target_drugs)} drugs in 'golden_set_library'...")
    
    # We fetch all and filter in python for case-insensitive matching or use 'in' if names match exactly
    # Let's try to fetch by name using 'in' with exact names first, if that fails we might need ilike loop
    
    res = supabase.table("golden_set_library")\
        .select("name, smiles_code, properties, relevance_score, ai_refined")\
        .in_("name", target_drugs)\
        .execute()
        
    found_map = {r['name']: r for r in res.data}
    
    # Write results to file to avoid encoding issues
    with open("db_results.txt", "w", encoding="utf-8") as f:
        f.write(f"{'Drug Name':<15} | {'SMILES':<10} | {'Target':<15} | {'Score':<5} | {'Refined':<7}\n")
        f.write("-" * 70 + "\n")
        
        for drug in target_drugs:
            record = found_map.get(drug)
            if not record:
                f.write(f"{drug:<15} | {'NOT FOUND':<10} | {'-':<15} | {'-':<5} | {'-':<7}\n")
                continue
                
            smiles_status = "✅ Found" if record.get('smiles_code') else "❌ Missing"
            
            props = record.get('properties', {})
            ai_analysis = props.get('ai_analysis', {})
            target = ai_analysis.get('target') or props.get('target') or "null"
            
            score = record.get('relevance_score', 0)
            refined = "YES" if record.get('ai_refined') else "NO"
            
            f.write(f"{drug:<15} | {smiles_status:<10} | {str(target)[:15]:<15} | {score:<5} | {refined:<7}\n")
            
            # Detailed check for Mylotarg
            if drug == "Mylotarg":
                f.write(f"\n   👉 Mylotarg Details:\n")
                f.write(f"      - Generic Name: {props.get('generic_name')}\n")
                f.write(f"      - SMILES: {record.get('smiles_code')[:50]}..." if record.get('smiles_code') else "      - SMILES: None\n")
                f.write(f"      - MoA Length: {len(props.get('mechanism_of_action', ''))}\n")
                f.write(f"      - Indication Length: {len(props.get('indication', ''))}\n")
    
    print("Verification complete. Results written to db_results.txt")

if __name__ == "__main__":
    asyncio.run(verify_db_results())
