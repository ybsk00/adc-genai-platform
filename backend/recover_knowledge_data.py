"""
Knowledge Base 불량 데이터 복구 스크립트
relevance_score = 0인 PubMed 데이터를 재분석하여 정상화

사용법:
    python recover_knowledge_data.py

주의: 실행 전 .env 환경 변수 설정 필요
"""
import os
import sys
import asyncio
import json
import logging
from datetime import datetime

# 프로젝트 루트 추가
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dotenv import load_dotenv
load_dotenv()

import google.generativeai as genai
from supabase import create_client
from json_repair import repair_json

# 로깅 설정
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Supabase 클라이언트
supabase = create_client(
    os.getenv("SUPABASE_URL"),
    os.getenv("SUPABASE_SERVICE_KEY")
)

# Gemini Safety Settings (의학 용어 차단 해제)
SAFETY_SETTINGS = [
    {"category": "HARM_CATEGORY_HATE_SPEECH", "threshold": "BLOCK_NONE"},
    {"category": "HARM_CATEGORY_HARASSMENT", "threshold": "BLOCK_NONE"},
    {"category": "HARM_CATEGORY_SEXUALLY_EXPLICIT", "threshold": "BLOCK_NONE"},
    {"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}
]


async def analyze_abstract(abstract: str, title: str) -> dict:
    """Gemini 2.0 Flash로 논문 분석"""
    system_prompt = """You are an expert ADC (Antibody-Drug Conjugate) researcher.
Analyze the provided scientific abstract and extract structured information.

Output MUST be a JSON object with these exact fields:
1. "target": Molecular target(s) mentioned (e.g., "HER2", "TROP2", "CD19"). Return null if not found.
2. "indication": Cancer type or disease indication. Return null if not found.
3. "summary": A concise 3-sentence summary focusing on clinical results.
4. "relevance_score": Float between 0.0 and 1.0 for ADC relevance.
5. "ai_reasoning": One sentence explaining the paper's importance for ADC research.

IMPORTANT: Return ONLY raw JSON. Do not use markdown formatting."""

    full_prompt = f"""{system_prompt}

Title: {title}
Abstract: {abstract[:3000]}"""

    try:
        genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))
        model = genai.GenerativeModel('gemini-2.5-flash')  # 2.5-flash (최신)
        
        response = model.generate_content(
            full_prompt,
            safety_settings=SAFETY_SETTINGS
        )
        
        content = response.text.strip()
        
        try:
            repaired = repair_json(content)
            return json.loads(repaired)
        except:
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0]
            elif "```" in content:
                content = content.split("```")[1].split("```")[0]
            return json.loads(content.strip())
            
    except Exception as e:
        logger.error(f"Gemini analysis error: {e}")
        return {
            "target": None,
            "indication": None,
            "summary": f"Analysis failed: {str(e)}",
            "relevance_score": 0.0,
            "ai_reasoning": f"Error: {str(e)}"
        }


async def recover_zero_score_data(batch_size: int = 50):
    """relevance_score = 0인 PubMed 데이터 복구"""
    
    logger.info("🔄 [Recovery] Starting knowledge_base data recovery...")
    
    # 1. 대상 데이터 조회
    result = supabase.table("knowledge_base")\
        .select("id, title, content")\
        .eq("source_type", "PubMed")\
        .eq("relevance_score", 0)\
        .limit(batch_size)\
        .execute()
    
    items = result.data
    
    if not items:
        logger.info("✨ No zero-score items to recover!")
        return {"recovered": 0}
    
    logger.info(f"📋 Found {len(items)} items to recover")
    
    recovered = 0
    errors = 0
    
    for idx, item in enumerate(items):
        try:
            logger.info(f"🔬 [{idx+1}/{len(items)}] Analyzing: {item['title'][:50]}...")
            
            # AI 분석
            analysis = await analyze_abstract(item["content"], item["title"])
            
            # 결과 검증
            if analysis.get("relevance_score", 0) == 0 and not analysis.get("target"):
                logger.warning(f"⚠️ Analysis still returned 0 for: {item['title'][:40]}")
            
            # DB 업데이트 (최소 필드만 - 존재하지 않는 컬럼 에러 방지)
            summary_text = f"Target: {analysis.get('target') or 'Unknown'}"
            if analysis.get("indication"):
                summary_text += f" | Indication: {analysis['indication']}"
            if analysis.get("summary"):
                summary_text += f"\n{analysis['summary']}"
            if analysis.get("ai_reasoning"):
                summary_text += f"\nReasoning: {analysis['ai_reasoning']}"
            
            # 최소 필드만 업데이트 (summary, relevance_score만 - rag_status, updated_at 제거)
            update_data = {
                "summary": summary_text[:1000],
                "relevance_score": float(analysis.get("relevance_score", 0.0))
            }
            
            supabase.table("knowledge_base").update(update_data).eq("id", item["id"]).execute()
            
            recovered += 1
            logger.info(f"✅ Recovered: Score={analysis.get('relevance_score', 0):.2f}, Target={analysis.get('target')}")
            
            # Rate limiting
            await asyncio.sleep(0.5)
            
        except Exception as e:
            errors += 1
            import traceback
            logger.error(f"❌ Recovery error for {item['id']}: {e}")
            logger.error(f"   Traceback: {traceback.format_exc()}")
    
    logger.info(f"🎉 Recovery complete! Recovered: {recovered}, Errors: {errors}")
    
    return {"recovered": recovered, "errors": errors}


async def verify_model():
    """Gemini 모델 연결 테스트"""
    logger.info("🔍 Testing Gemini connection...")
    
    try:
        genai.configure(api_key=os.getenv("GOOGLE_API_KEY"))
        model = genai.GenerativeModel('gemini-2.5-flash')  # 2.5-flash
        
        response = model.generate_content(
            "Say 'OK' if you can respond.",
            safety_settings=SAFETY_SETTINGS
        )
        
        if response.text:
            logger.info(f"✅ Gemini connection OK: {response.text.strip()}")
            return True
        else:
            logger.error("❌ Empty response from Gemini")
            return False
            
    except Exception as e:
        logger.error(f"❌ Gemini connection failed: {e}")
        return False


async def main():
    import argparse
    parser = argparse.ArgumentParser(description="Knowledge Base Data Recovery")
    parser.add_argument("--batch", type=int, default=50, help="Batch size (default: 50)")
    parser.add_argument("--test", action="store_true", help="Test Gemini connection only")
    args = parser.parse_args()
    
    if args.test:
        await verify_model()
    else:
        # 1. 먼저 모델 테스트
        if not await verify_model():
            logger.error("❌ Cannot proceed - Gemini connection failed")
            return
        
        # 2. 데이터 복구 실행
        result = await recover_zero_score_data(args.batch)
        print(f"\n📊 Result: {result}")


if __name__ == "__main__":
    asyncio.run(main())
