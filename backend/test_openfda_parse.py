"""
OpenFDA AI Refiner 단일 레코드 테스트
ELAHERE, ZYNLONTA에 대해 텍스트 파싱 및 AI 분석 테스트
"""
import asyncio
import sys
sys.path.insert(0, ".")

from app.core.supabase import supabase
from app.services.ai_refiner import ai_refiner
import json

async def test_single_record(drug_name: str):
    print(f"\n{'='*60}")
    print(f"🔍 Testing: {drug_name}")
    print('='*60)
    
    # 1. DB에서 레코드 조회
    result = supabase.table("golden_set_library")\
        .select("*")\
        .ilike("name", f"%{drug_name}%")\
        .eq("enrichment_source", "open_fda_api")\
        .limit(1)\
        .execute()
    
    if not result.data:
        print(f"❌ No record found for {drug_name}")
        return
    
    record = result.data[0]
    print(f"✅ Found record ID: {record['id']}")
    print(f"   Name: {record['name']}")
    
    # 2. Properties 구조 확인
    properties = record.get("properties", {})
    print(f"\n📂 Properties Keys: {list(properties.keys())}")
    
    if "fda_label" in properties:
        fda_label = properties["fda_label"]
        print(f"   fda_label Keys: {list(fda_label.keys())}")
        
        # 주요 필드 출력
        for key in ["indications_and_usage", "indication", "mechanism_of_action", "generic_name"]:
            if key in fda_label:
                val = fda_label[key]
                print(f"   - {key}: {val[:100] if val else 'EMPTY'}...")
    else:
        # fda_label 없을 때 properties 직접 확인
        for key in ["indications_and_usage", "indication", "mechanism_of_action", "description"]:
            if key in properties:
                val = properties[key]
                print(f"   - {key}: {val[:100] if val else 'EMPTY'}...")
    
    # 3. AI Refiner 실행
    print(f"\n🤖 Running AI Refiner...")
    analysis = await ai_refiner.refine_single_record(record)
    
    if analysis:
        print(f"\n📊 AI Analysis Result:")
        print(json.dumps(analysis, indent=2, ensure_ascii=False))
    else:
        print(f"❌ Analysis failed!")
    
    return analysis

async def main():
    # 테스트 대상 약물
    drugs = ["ELAHERE", "ZYNLONTA", "Mylotarg"]
    
    for drug in drugs:
        await test_single_record(drug)
        print()
    
    print("\n" + "="*60)
    print("✅ Test Complete!")

if __name__ == "__main__":
    asyncio.run(main())
