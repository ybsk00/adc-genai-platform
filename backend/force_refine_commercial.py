
import asyncio
import logging
from app.core.supabase import supabase
from app.services.ai_refiner import ai_refiner

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def force_refine_commercial():
    print("🚀 Force Refinement Start: Commercial Reagents")
    
    # 1. 미정제 데이터 조회
    res = supabase.table("commercial_reagents")\
        .select("*")\
        .eq("ai_refined", False)\
        .limit(20)\
        .execute()
    
    if not res.data:
        print("✨ No pending records in commercial_reagents.")
        return

    print(f"🔍 Found {len(res.data)} pending records. Starting AI...")
    
    for record in res.data:
        print(f"🤖 Processing: {record.get('product_name')}")
        analysis = await ai_refiner.refine_single_record(record)
        
        if analysis and "error" not in analysis:
            update_data = {
                "target": analysis.get("target"),
                "ai_refined": True,
                "properties": {**record.get("properties", {}), "ai_analysis": analysis}
            }
            supabase.table("commercial_reagents").update(update_data).eq("id", record["id"]).execute()
            print(f"✅ Success: {record.get('product_name')}")
        else:
            print(f"❌ Failed: {record.get('product_name')} - {analysis.get('error') if analysis else 'Unknown'}")
        
        await asyncio.sleep(1)

if __name__ == "__main__":
    asyncio.run(force_refine_commercial())
