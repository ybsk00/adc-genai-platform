"""
PubMed/BioRxiv 데이터 AI 정제 서비스
Gemini 2.0 Flash를 사용하여 초록을 분석하고 요약, 관련성 점수, 근거를 생성함.
"""
import json
import logging
import asyncio
from typing import List, Dict, Any, Optional
from datetime import datetime

from app.core.config import settings
from app.core.supabase import supabase
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.prompts import ChatPromptTemplate

logger = logging.getLogger(__name__)

class KnowledgeBaseRefiner:
    def __init__(self):
        # Gemini 2.0 Flash (빠르고 저렴함)
        self.llm = ChatGoogleGenerativeAI(
            model="gemini-2.0-flash-exp", # 또는 사용 가능한 최신 Flash 모델
            google_api_key=settings.GOOGLE_API_KEY,
            temperature=0,
            convert_system_message_to_human=True
        )

    async def process_pending_items(self, job_id: Optional[str] = None, batch_size: int = 20):
        """rag_status='pending'인 항목을 가져와서 AI 분석 수행"""
        try:
            # 1. Fetch pending items
            items_res = supabase.table("knowledge_base")\
                .select("*")\
                .eq("rag_status", "pending")\
                .limit(batch_size)\
                .execute()
            
            items = items_res.data
            if not items:
                logger.info("No pending items to refine.")
                return 0

            processed_count = 0
            
            # 2. Process each item
            for item in items:
                try:
                    result = await self.analyze_abstract(item["content"])
                    
                    # 3. Update DB
                    supabase.table("knowledge_base").update({
                        "summary": result.get("summary", ""),
                        "relevance_score": result.get("relevance_score", 0.0),
                        "ai_reasoning": result.get("ai_reasoning", ""),
                        "rag_status": "processed",
                        "updated_at": datetime.utcnow().isoformat()
                    }).eq("id", item["id"]).execute()
                    
                    processed_count += 1
                    logger.info(f"✅ Refined Knowledge Base Item: {item['id']}")
                    
                    # Rate limiting (Gemini Flash is fast but still good to be safe)
                    await asyncio.sleep(0.5)
                    
                except Exception as e:
                    logger.error(f"Failed to refine item {item['id']}: {e}")
                    # 실패 시 status를 error로 바꾸거나 retry count를 늘리는 로직이 있으면 좋음
                    # 여기서는 일단 pending으로 유지 (다음 실행 때 재시도)
            
            return processed_count

        except Exception as e:
            logger.error(f"Knowledge Refiner Error: {e}")
            return 0

    async def analyze_abstract(self, abstract_text: str) -> Dict[str, Any]:
        """Gemini Flash를 사용하여 초록 분석 (JSON 출력 강제)"""
        
        system_prompt = """You are an expert researcher in Antibody-Drug Conjugates (ADC) for oncology.
Analyze the provided scientific abstract.

Output MUST be a JSON object with these fields:
1. "summary": A concise 3-sentence summary focusing on the target antigen, payload, and linker technology.
2. "relevance_score": A float between 0.0 and 1.0.
   - 1.0: Directly about ADC drug development or clinical trials for cancer.
   - 0.5: General antibody research or non-ADC conjugates.
   - 0.1: Irrelevant (e.g., diagnostic imaging, viral vaccines).
3. "ai_reasoning": A short explanation of why you gave this score.

IMPORTANT: Return ONLY raw JSON. Do not use markdown formatting like ```json ... ```.
"""
        
        user_prompt = f"Abstract:\n{abstract_text}"
        
        try:
            # LangChain의 bind를 사용하여 response_mime_type 설정 (Gemini 전용)
            # 또는 프롬프트 엔지니어링으로 해결. 
            # langchain-google-genai 최신 버전에서는 bind(generation_config=...) 지원
            
            prompt = ChatPromptTemplate.from_messages([
                ("system", system_prompt),
                ("user", user_prompt)
            ])
            
            chain = prompt | self.llm
            
            # JSON 모드 강제 (모델 파라미터 지원 시)
            # 현재 라이브러리 버전에 따라 다를 수 있으므로 프롬프트에 의존하되, 
            # 파싱 로직을 견고하게 작성
            
            response = await chain.ainvoke({})
            content = response.content.strip()
            
            # Markdown code block 제거
            if content.startswith("```json"):
                content = content[7:]
            if content.endswith("```"):
                content = content[:-3]
            
            return json.loads(content.strip())
            
        except Exception as e:
            logger.error(f"AI Analysis Error: {e}")
            return {
                "summary": "Analysis Failed",
                "relevance_score": 0.0,
                "ai_reasoning": f"Error: {str(e)}"
            }

knowledge_refiner = KnowledgeBaseRefiner()
