import asyncio
import logging
import httpx
from app.services.openfda_service import openfda_service
from app.core.supabase import supabase

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("FORCE_REFETCH")

async def force_refetch():
    print("\n" + "="*80)
    print("🚀 [FORCE REFETCH] Updating 5 Key ADC Drugs from OpenFDA")
    print("="*80 + "\n")

    targets = [
        {"name": "Mylotarg", "query": 'openfda.brand_name:"MYLOTARG"'},
        {"name": "ELAHERE", "query": 'openfda.brand_name:"ELAHERE"'},
        {"name": "ZYNLONTA", "query": 'openfda.brand_name:"ZYNLONTA"'},
        {"name": "Besponsa", "query": 'openfda.brand_name:"BESPONSA"'},
        {"name": "Blenrep", "query": 'openfda.brand_name:"BLENREP"'}
    ]

    async with httpx.AsyncClient(timeout=30.0) as client:
        for target in targets:
            name = target["name"]
            query = target["query"]
            print(f"🔄 Fetching {name} (Query: {query})...")
            
            # 1. Fetch from OpenFDA
            url = f"https://api.fda.gov/drug/label.json?search={query}&limit=1"
            try:
                res = await client.get(url)
                if res.status_code != 200:
                    print(f"   ❌ Failed to fetch {name}: {res.status_code}")
                    continue
                
                data = res.json()
                results = data.get("results", [])
                if not results:
                    print(f"   ⚠️ No results found for {name}")
                    continue
                
                label_data = results[0]
                
                # 2. Extract Info (using the FIXED extract_golden_info)
                golden_data = openfda_service.extract_golden_info(label_data)
                
                # 3. Verify Text Lengths immediately
                props = golden_data["properties"]
                ind_len = len(props.get("indication", ""))
                moa_len = len(props.get("mechanism_of_action", ""))
                
                print(f"   ✅ Fetched Data:")
                print(f"      - Indication Length: {ind_len}")
                print(f"      - MoA Length: {moa_len}")
                
                if ind_len < 100 or moa_len < 100:
                    print("      ⚠️ WARNING: Text seems short!")
                
                # 4. Update DB (Force Update)
                # We want to ensure this data overwrites/merges correctly.
                # We will fetch the existing ID first to do an update.
                
                existing = supabase.table("golden_set_library").select("id, properties").ilike("name", name).execute()
                
                if existing.data:
                    rec_id = existing.data[0]["id"]
                    existing_props = existing.data[0].get("properties", {}) or {}
                    
                    # Merge logic: Put new data into fda_label, but ALSO update top-level keys if they are empty/short?
                    # Actually, let's just trust the fda_label approach since we fixed ai_refiner to look there.
                    # BUT, to be absolutely sure as per user request, let's put the text in fda_label AND ensure it's accessible.
                    
                    merged_props = existing_props.copy()
                    merged_props["fda_label"] = props # The whole properties from extract_golden_info IS the fda_label content basically
                    
                    # Also update top-level keys for visibility if needed, but let's stick to the pattern
                    # extract_golden_info returns a dict where 'properties' contains 'indication', 'mechanism_of_action' etc.
                    
                    supabase.table("golden_set_library").update({
                        "properties": merged_props,
                        "ai_refined": False, # Reset so AI runs again
                        "relevance_score": 0, # Reset score
                        "updated_at": "now()"
                    }).eq("id", rec_id).execute()
                    
                    print(f"   💾 Updated DB Record: {rec_id}")
                    
                else:
                    # Insert new
                    print(f"   ⚠️ Record not found in DB, inserting new...")
                    supabase.table("golden_set_library").insert(golden_data).execute()
                    print(f"   💾 Inserted New Record")

            except Exception as e:
                print(f"   ❌ Error processing {name}: {e}")

    print("\n" + "="*80)
    print("✅ Refetch Complete. Please run AI Refiner (Run Full) now.")

if __name__ == "__main__":
    asyncio.run(force_refetch())
