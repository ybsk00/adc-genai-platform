"""
OpenFDA AI Refiner 하드코딩 테스트 (강제 주입 테스트)
DB나 파싱 로직을 거치지 않고, 텍스트를 직접 주입하여 AI 로직만 검증합니다.
"""
import asyncio
import sys
import json
import logging

# 로깅 설정
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

sys.path.insert(0, ".")

from app.services.ai_refiner import ai_refiner

async def test_hardcoded_elahere():
    print("\n" + "="*60)
    print("🧪 ELAHERE 하드코딩 테스트 (Indication 텍스트 강제 주입)")
    print("="*60)
    
    # ELAHERE의 실제 Indication 텍스트 (OpenFDA)
    indication_text = (
        "ELAHERE is indicated for the treatment of adult patients with folate receptor-alpha (FRα) positive, "
        "platinum-resistant epithelial ovarian, fallopian tube, or primary peritoneal cancer, "
        "who have received one to three prior systemic treatment regimens. "
        "Select patients for therapy based on an FDA-approved test."
    )
    
    # 가짜 레코드 생성 (MoA는 일부러 비움)
    fake_record = {
        "name": "ELAHERE",
        "enrichment_source": "open_fda_api",
        "properties": {
            "fda_label": {
                "indication": indication_text,
                "mechanism_of_action": "", # MoA 비움 -> Indication에서 찾아야 함
                "generic_name": "MIRVETUXIMAB SORAVTANSINE",
                "boxed_warning": "WARNING: OCULAR TOXICITY"
            }
        }
    }
    
    print(f"📝 Input Data:")
    print(f"   Name: {fake_record['name']}")
    print(f"   Indication: {indication_text[:100]}...")
    print(f"   MoA: (Empty)")
    print("-" * 60)
    
    # AI Refiner 실행
    print("🤖 Running AI Refiner...")
    try:
        analysis = await ai_refiner.refine_single_record(fake_record)
        
        print("\n📊 AI Analysis Result:")
        print(json.dumps(analysis, indent=2, ensure_ascii=False))
        
        # 검증
        target = analysis.get("target", "")
        score = analysis.get("relevance_score", 0)
        
        print("\n" + "="*60)
        print("🧐 검증 결과:")
        
        if "FRα" in target or "folate receptor-alpha" in target:
            print("✅ Target Extraction: SUCCESS (Found FRα)")
        else:
            print(f"❌ Target Extraction: FAILED (Got '{target}')")
            
        if score >= 0.8:
             print(f"✅ Relevance Score: SUCCESS ({score})")
        else:
             print(f"❌ Relevance Score: FAILED ({score})")
             
    except Exception as e:
        print(f"❌ Error during execution: {e}")

if __name__ == "__main__":
    asyncio.run(test_hardcoded_elahere())
