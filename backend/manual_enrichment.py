import asyncio
import argparse
import json
import os
import sys

# Add backend directory to path to import app modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from app.core.supabase import supabase

async def manual_update(target_name: str, updates: dict):
    """
    Manually update a record in golden_set_library.
    """
    print(f"🔍 Searching for: {target_name}")
    
    # 1. Find the record (Case-insensitive search)
    res = supabase.table("golden_set_library")\
        .select("id, name, properties")\
        .ilike("name", f"%{target_name}%")\
        .limit(5)\
        .execute()
        
    if not res.data:
        print(f"❌ No record found matching '{target_name}'")
        return

    if len(res.data) > 1:
        print(f"⚠️ Multiple records found ({len(res.data)}). Please be more specific.")
        for item in res.data:
            print(f" - {item['name']} (ID: {item['id']})")
        return

    record = res.data[0]
    record_id = record['id']
    print(f"✅ Found record: {record['name']} (ID: {record_id})")
    
    # 2. Prepare update payload
    payload = {
        "enrichment_source": "Manual Review", # Mark as manually reviewed
        "confidence_score": 1.0,              # Manual = 100% confidence
        "ai_refined": True,
        "review_required": False              # Clear review flag
    }
    
    # Top-level columns in golden_set_library
    TOP_LEVEL_FIELDS = [
        "outcome_type", "dar", "orr_pct", "os_months", "pfs_months", 
        "dor_months", "patient_count", "adverse_events_grade3_pct",
        "target_1", "target_symbol", "antibody_format", "payload_smiles"
    ]
    
    properties = record.get("properties", {}) or {}
    
    print(f"📝 Applying updates: {updates}")
    
    for key, value in updates.items():
        if key in TOP_LEVEL_FIELDS:
            payload[key] = value
        else:
            # Update properties for non-column fields
            properties[key] = value
            
    payload["properties"] = properties
    
    # 3. Execute Update
    try:
        print("⏳ Sending update request to Supabase...")
        data = supabase.table("golden_set_library").update(payload).eq("id", record_id).execute()
        
        if not data.data:
            print("❌ Update failed: No data returned from database.")
            sys.exit(1)
            
        updated_record = data.data[0]
        print("✅ Database confirmed update.")
        
        # 4. Strict Verification
        print("\n� Verifying updated fields...")
        failed_checks = []
        
        for key, expected_value in updates.items():
            # Handle nested properties vs top-level columns
            actual_value = updated_record.get(key)
            if key not in TOP_LEVEL_FIELDS:
                actual_value = updated_record.get("properties", {}).get(key)
            
            # Simple equality check (convert to str for comparison to avoid type mismatches like float vs int)
            if str(actual_value) != str(expected_value):
                failed_checks.append(f"   - {key}: Expected '{expected_value}', Got '{actual_value}'")
        
        if failed_checks:
            print("❌ VERIFICATION FAILED: The following fields were not updated correctly:")
            for failure in failed_checks:
                print(failure)
            print("🛑 Stopping execution due to verification failure.")
            sys.exit(1)
        else:
            print("✨ SUCCESS: All fields verified correctly!")
            print(json.dumps(updated_record, indent=2, default=str))
            
    except Exception as e:
        print(f"❌ Update failed with exception: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Manual Enrichment Tool for ADC Platform")
    parser.add_argument("--name", required=True, help="Name of the drug/trial to update")
    
    # Optional fields to update
    parser.add_argument("--orr", type=float, help="ORR %")
    parser.add_argument("--os", type=float, help="OS (months)")
    parser.add_argument("--pfs", type=float, help="PFS (months)")
    parser.add_argument("--dar", type=float, help="Drug-to-Antibody Ratio")
    parser.add_argument("--ae", type=float, help="Grade 3+ AE %")
    parser.add_argument("--count", type=int, help="Patient Count")
    parser.add_argument("--target", type=str, help="Target Symbol (e.g. ERBB2)")
    parser.add_argument("--smiles", type=str, help="Payload SMILES")
    parser.add_argument("--outcome", type=str, help="Outcome Type (Success/Failure)")
    
    # Extended fields
    parser.add_argument("--format", type=str, help="Antibody Format (e.g. IgG1)")
    parser.add_argument("--linker", type=str, help="Linker Type")
    parser.add_argument("--gene_id", type=str, help="NCBI Gene ID")
    parser.add_argument("--host", type=str, help="Host Species")
    parser.add_argument("--uniprot", type=str, help="UniProt ID")
    parser.add_argument("--affinity", type=str, help="Binding Affinity")
    parser.add_argument("--isotype", type=str, help="Isotype")
    parser.add_argument("--mw", type=str, help="Molecular Weight")
    
    args = parser.parse_args()
    
    # Construct updates dictionary
    updates = {}
    if args.orr is not None: updates["orr_pct"] = args.orr
    if args.os is not None: updates["os_months"] = args.os
    if args.pfs is not None: updates["pfs_months"] = args.pfs
    if args.dar is not None: updates["dar"] = args.dar
    if args.ae is not None: updates["adverse_events_grade3_pct"] = args.ae
    if args.count is not None: updates["patient_count"] = args.count
    if args.target is not None: updates["target_symbol"] = args.target
    if args.smiles is not None: updates["payload_smiles"] = args.smiles
    if args.outcome is not None: updates["outcome_type"] = args.outcome
    
    if args.format is not None: updates["antibody_format"] = args.format
    if args.linker is not None: updates["linker_type"] = args.linker
    if args.gene_id is not None: updates["gene_id"] = args.gene_id
    if args.host is not None: updates["host_species"] = args.host
    if args.uniprot is not None: updates["uniprot_id"] = args.uniprot
    if args.affinity is not None: updates["binding_affinity"] = args.affinity
    if args.isotype is not None: updates["isotype"] = args.isotype
    if args.mw is not None: updates["molecular_weight"] = args.mw
    
    if not updates:
        print("⚠️ No updates specified. Use flags like --orr, --os, --dar etc.")
        return

    # Run async function
    asyncio.run(manual_update(args.name, updates))

if __name__ == "__main__":
    main()
