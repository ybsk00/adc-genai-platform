import asyncio
import json
import logging
import sys
import os
import traceback
from dotenv import load_dotenv

from app.core.supabase import supabase
from app.services.ai_refiner import ai_refiner
from app.services.knowledge_refiner import knowledge_refiner

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

async def test_clinical_refiner_single():
    print("\n=== [1] Clinical Trials Refiner Test (AIRefiner) ===")
    # 전체 카운트 확인
    count_res = supabase.table("golden_set_library").select("id", count="exact").execute()
    print(f"Total records in golden_set_library: {count_res.count}")

    # 모든 레코드의 에러 상태 확인 (디버깅용)
    all_res = supabase.table("golden_set_library").select("id, processing_error, ai_refined").limit(10).execute()
    print(f"Sample Records: {all_res.data}")

    # 실패한 레코드 하나 가져오기 (단순화된 쿼리)
    res = supabase.table("golden_set_library")\
        .select("*")\
        .not_.is_("processing_error", "null")\
        .limit(1)\
        .execute()
    
    if not res.data:
        print("실패한 Clinical Trials 레코드가 없습니다. (이미 리셋되었거나 아직 처리 전일 수 있음)")
        # 대기 중인 레코드라도 하나 가져오기
        res = supabase.table("golden_set_library").select("*").eq("ai_refined", False).limit(1).execute()
        if not res.data:
            print("테스트할 레코드가 아예 없습니다.")
            return

    record = res.data[0]
    print(f"Target Record ID: {record['id']}")
    print(f"Target Name: {record.get('name')}")
    
    # ai_refiner.refine_single_record의 로직을 직접 실행하며 로그 확인
    try:
        # AIRefiner의 내부 로직을 모방하되 로그를 더 자세히 찍음
        from langchain_openai import ChatOpenAI
        from app.core.config import settings
        from langchain_core.prompts import ChatPromptTemplate
        
        llm = ChatOpenAI(
            model=settings.FAST_LLM,
            temperature=0,
            api_key=settings.OPENAI_API_KEY
        )
        
        properties = record.get("properties", {})
        title = record.get("name", "No Title")
        description = properties.get("brief_summary", "")
        overall_status = properties.get("overall_status", "")
        why_stopped = properties.get("why_stopped", "")
        phase = properties.get("phase", "")
        
        system_prompt = """You are a Clinical Trial Analyst specializing in ADC research.
Output ONLY valid JSON in this exact format:
{
    "drug_name": "extracted drug/compound name or null",
    "target": "molecular target (e.g., HER2, TROP2) or null",
    "outcome_type": "Success|Failure|Ongoing|Unknown",
    "failure_reason": "reason if failed, null otherwise",
    "confidence": 0.0-1.0
}"""
        user_prompt = f"Title: {title}\nPhase: {phase}\nStatus: {overall_status}\nWhy Stopped: {why_stopped}\nDescription: {description[:1000]}"
        
        prompt = ChatPromptTemplate.from_messages([("system", system_prompt), ("user", user_prompt)])
        chain = prompt | llm
        
        print("\n--- LLM Requesting... ---")
        response = await chain.ainvoke({})
        content = response.content.strip()
        print(f"Raw LLM Response:\n{content}")
        
        print("\n--- Parsing Attempt... ---")
        # 기존 파싱 로직 적용
        parsed_content = content
        if "```json" in content:
            parsed_content = content.split("```json")[1].split("```")[0]
        elif "```" in content:
            parsed_content = content.split("```")[1].split("```")[0]
        
        print(f"Cleaned Content for JSON Load:\n{parsed_content}")
        
        try:
            analysis = json.loads(parsed_content.strip())
            print("✅ JSON Parsing Success!")
            print(json.dumps(analysis, indent=2))
        except Exception as e:
            print(f"❌ JSON Parsing Failed: {e}")
            
    except Exception as e:
        print(f"❌ Test Error: {e}")

async def test_knowledge_refiner_single():
    print("\n=== [2] Knowledge Base Refiner Test (Gemini) ===")
    # 실패한 레코드 하나 가져오기 (summary가 'Analysis Failed'인 것)
    res = supabase.table("knowledge_base")\
        .select("*")\
        .eq("summary", "Analysis Failed")\
        .limit(1)\
        .execute()
    
    if not res.data:
        print("실패한 Knowledge Base 레코드가 없습니다.")
        # 대기 중인 레코드라도 하나 가져오기
        res = supabase.table("knowledge_base").select("*").eq("rag_status", "pending").limit(1).execute()
        if not res.data:
            print("테스트할 레코드가 아예 없습니다.")
            return

    item = res.data[0]
    print(f"Target Item ID: {item['id']}")
    print(f"Target Title: {item.get('title')}")
    
    try:
        # KnowledgeBaseRefiner의 analyze_abstract 로직 실행
        result = await knowledge_refiner.analyze_abstract(item["content"])
        print("\n--- Refiner Result ---")
        print(json.dumps(result, indent=2))
        
        if result.get("summary") == "Analysis Failed":
            print("❌ 여전히 실패함. 원인:", result.get("ai_reasoning"))
        else:
            print("✅ 분석 성공!")
            
    except Exception as e:
        print(f"❌ Test Error: {e}")

if __name__ == "__main__":
    asyncio.run(test_clinical_refiner_single())
    asyncio.run(test_knowledge_refiner_single())
