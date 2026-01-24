import asyncio
from app.core.supabase import supabase

async def verify_db():
    print("🔍 Verifying Ambeed Data in DB...")
    
    # Query for source_name = 'Ambeed'
    res = supabase.table('commercial_reagents').select('*').eq('source_name', 'Ambeed').limit(5).execute()
    
    if res.data:
        print(f"✅ Found {len(res.data)} records from Ambeed.")
        for item in res.data:
            print(f"   - {item.get('product_name')} (Cat: {item.get('ambeed_cat_no')})")
            print(f"     Target: {item.get('target')}")
            print(f"     SMILES: {item.get('smiles_code')}")
    else:
        print("❌ No records found for Ambeed.")

if __name__ == "__main__":
    asyncio.run(verify_db())
