import os
import sys
import asyncio
import logging
import json
from datetime import datetime

# 1. Load .env manually to ensure DB connection and API keys work
env_path = os.path.join(os.getcwd(), "backend", ".env")
if os.path.exists(env_path):
    with open(env_path, "r", encoding="utf-8") as f:
        for line in f:
            if "=" in line and not line.strip().startswith("#"):
                key, value = line.strip().split("=", 1)
                value = value.strip().strip("'").strip('"')
                os.environ[key] = value
else:
    print("⚠️ Warning: .env not found")

# 2. Add backend to sys.path
sys.path.append(os.path.join(os.getcwd(), "backend"))

from app.core.supabase import supabase
from app.services.ai_refiner import ai_refiner

# Configure Logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DEBUG_PIPELINE")

async def test_full_pipeline():
    print("🚀 Starting Debug Test for Full Pipeline (DB Save + AI Trigger)...")
    
    # 3. Create dummy data for AI Test
    # Using a real-looking ADC payload name to check if AI identifies it correctly
    test_id = "DEBUG-AI-TEST-002"
    test_data = {
        "ambeed_cat_no": test_id,
        "cas_number": "123-45-6", # Fake CAS
        "product_name": "MMAE (Monomethyl auristatin E)", # Known ADC payload
        "product_url": "https://example.com/mmae-debug",
        "category": "ADC Toxin",
        "smiles_code": None, # Expect AI/PubChem to try to fill this, or at least analyze
        "target": None, # Expect AI to fill
        "properties": {"test": "pipeline_check"},
        "summary": "MMAE is a potent antimitotic agent used as a payload in antibody-drug conjugates (ADCs). It targets Tubulin.",
        "price_data": [],
        "source_name": "Creative Biolabs",
        "crawled_at": datetime.utcnow().isoformat()
    }

    print(f"   📝 Step 1: Inserting Test Record ({test_id})...")
    
    try:
        # Save to DB (The part we fixed)
        res = supabase.table("commercial_reagents").upsert(test_data, on_conflict="ambeed_cat_no").execute()
        
        if not res.data:
            print("   ❌ DB Save Failed (No data returned). Pipeline aborted.")
            return

        record = res.data[0]
        record_id = record['id']
        print(f"   ✅ DB Save Successful! ID: {record_id}")
        
        # 4. Trigger AI Refiner Manually (Simulating the async call in crawler)
        print(f"\n   🧠 Step 2: Triggering AI Refiner for ID: {record_id}...")
        print("      (This calls Gemini, so it might take a few seconds)")
        
        # Note: refine_single_record returns the analysis dict but doesn't update DB itself 
        # (The crawler calls update after getting the result).
        # Wait, let's check creative_biolabs_crawler.py again.
        # It calls: analysis = await ai_refiner.refine_single_record(record)
        # Then: supabase.table("commercial_reagents").update(...).eq("id", record["id"]).execute()
        
        # So I will replicate that logic here to prove it works.
        
        start_time = datetime.now()
        analysis = await ai_refiner.refine_single_record(record)
        duration = (datetime.now() - start_time).total_seconds()
        
        print(f"   ⏱️ Analysis finished in {duration:.2f}s")
        
        if analysis and "error" not in analysis:
            print("\n   ✨ AI Analysis Result:")
            print(json.dumps(analysis, indent=2))
            
            # Verify critical fields
            target = analysis.get("target")
            category = analysis.get("category")
            
            if target or category:
                print(f"      ✅ Extracted Target: {target}")
                print(f"      ✅ Extracted Category: {category}")
            else:
                print("      ⚠️ AI returned result but Target/Category might be empty.")
                
            # Simulate the DB Update (Final Step)
            print("\n   💾 Step 3: Updating DB with Analysis...")
            update_data = {
                "target": analysis.get("target"),
                "ai_refined": True,
                "properties": {**record.get("properties", {}), "ai_analysis": analysis}
            }
            
            upd_res = supabase.table("commercial_reagents").update(update_data).eq("id", record_id).execute()
            
            # Verify the final state
            final_check = supabase.table("commercial_reagents").select("*").eq("id", record_id).execute()
            final_record = final_check.data[0]
            
            if final_record.get('ai_refined') == True:
                print(f"   ✅ Final DB Verification: ai_refined = True")
                print("   🎉 Full Pipeline Test SUCCESS!")
            else:
                print("   ❌ Final DB Verification Failed: ai_refined is still False")

        else:
            print(f"   ❌ AI Analysis Failed or Returned Error: {analysis}")

    except Exception as e:
        print(f"\n   ❌ Test Failed with Exception: {e}")

if __name__ == "__main__":
    asyncio.run(test_full_pipeline())
