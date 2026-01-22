"""
ELAHERE 최종 프롬프트 전문 추출기
Gemini에게 전달되는 실제 텍스트를 파일로 저장
"""
import sys
sys.path.insert(0, ".")

from app.core.supabase import supabase
import json

def extract_prompt_for_elahere():
    print("="*60)
    print("🔍 ELAHERE 최종 프롬프트 추출")
    print("="*60)
    
    # 1. DB에서 ELAHERE 레코드 조회
    result = supabase.table("golden_set_library")\
        .select("*")\
        .ilike("name", "%ELAHERE%")\
        .eq("enrichment_source", "open_fda_api")\
        .limit(1)\
        .execute()
    
    if not result.data:
        print("❌ ELAHERE not found!")
        return
    
    record = result.data[0]
    properties = record.get("properties", {})
    title = record.get("name", "No Title")
    
    print(f"\n📂 Record ID: {record['id']}")
    print(f"📂 Name: {title}")
    
    # 2. Properties 구조 상세 출력
    print(f"\n{'='*40}")
    print("📋 PROPERTIES 구조 (전체)")
    print("="*40)
    print(json.dumps(properties, indent=2, ensure_ascii=False)[:5000])
    print("\n... (truncated if too long)")
    
    # 3. fda_label 확인
    fda_label = properties.get("fda_label", {})
    print(f"\n{'='*40}")
    print("📋 fda_label 키 목록:")
    print("="*40)
    for key in fda_label.keys():
        val = fda_label[key]
        val_type = type(val).__name__
        if isinstance(val, list):
            print(f"  - {key}: LIST ({len(val)} items)")
            if val:
                print(f"       첫 번째 항목: {str(val[0])[:100]}...")
        elif isinstance(val, str):
            print(f"  - {key}: STRING ({len(val)} chars)")
            print(f"       내용: {val[:100]}...")
        else:
            print(f"  - {key}: {val_type}")
    
    # 4. 텍스트 추출 시도 (리스트 처리 포함)
    print(f"\n{'='*40}")
    print("🔧 텍스트 추출 시도 (리스트 → 문자열 변환)")
    print("="*40)
    
    def safe_extract(data, key):
        """리스트면 join, 문자열이면 그대로"""
        val = data.get(key, "")
        if isinstance(val, list):
            return " ".join(str(v) for v in val)
        return val or ""
    
    description = (
        safe_extract(fda_label, "indications_and_usage") or
        safe_extract(fda_label, "indication") or
        safe_extract(properties, "indications_and_usage") or
        safe_extract(properties, "indication") or
        safe_extract(fda_label, "description") or
        safe_extract(properties, "description")
    )
    
    moa = safe_extract(fda_label, "mechanism_of_action") or safe_extract(properties, "mechanism_of_action")
    boxed_warning = safe_extract(fda_label, "boxed_warning") or safe_extract(fda_label, "warnings")
    generic_name = safe_extract(fda_label, "generic_name") or safe_extract(properties, "generic_name")
    
    print(f"Description 추출 ({len(description)} chars):")
    print(f"  → {description[:300]}..." if description else "  → EMPTY!")
    print(f"\nMoA 추출 ({len(moa)} chars):")
    print(f"  → {moa[:300]}..." if moa else "  → EMPTY!")
    print(f"\nGeneric Name: {generic_name or 'EMPTY!'}")
    
    # 5. 최종 프롬프트 생성
    system_prompt = """You are a Pharmaceutical Regulatory Affairs Specialist.
Analyze the FDA Drug Label data for an ADC (Antibody-Drug Conjugate).

Output ONLY valid JSON:
{
    "drug_name": "extracted drug name",
    "target": "molecular target (e.g., HER2, CD19, FRα, TROP2) or null",
    "outcome_type": "Success",
    "approval_status": "Approved",
    "boxed_warning": "Summary of Boxed Warning or 'None'",
    "indication": "Primary indication (e.g., Breast Cancer)",
    "relevance_score": 1.0,
    "confidence": 0.0-1.0
}
"""
    full_prompt = f"""{system_prompt}

FDA Label Data:
Name: {title}
Generic Name: {generic_name}
Indication: {description[:500] if description else "N/A"}
Mechanism of Action: {moa[:800] if moa else "N/A"}
Boxed Warning: {boxed_warning[:300] if boxed_warning else "N/A"}
"""
    
    # 6. 최종 프롬프트 파일 저장
    output_file = "ELAHERE_final_prompt.txt"
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("="*60 + "\n")
        f.write("ELAHERE - Gemini 최종 프롬프트 전문\n")
        f.write("="*60 + "\n\n")
        f.write(full_prompt)
        f.write("\n\n" + "="*60 + "\n")
        f.write(f"프롬프트 총 길이: {len(full_prompt)} chars\n")
        f.write("="*60 + "\n")
    
    print(f"\n{'='*40}")
    print(f"✅ 최종 프롬프트 저장 완료: {output_file}")
    print(f"   프롬프트 총 길이: {len(full_prompt)} chars")
    print("="*40)
    
    # 7. 콘솔에도 전체 출력
    print(f"\n{'='*60}")
    print("📤 최종 프롬프트 전문 (Gemini에게 전달되는 내용)")
    print("="*60)
    print(full_prompt)

if __name__ == "__main__":
    extract_prompt_for_elahere()
