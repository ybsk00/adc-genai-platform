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
import google.generativeai as genai
from app.services.cost_tracker import cost_tracker
from json_repair import repair_json

logger = logging.getLogger(__name__)

class KnowledgeBaseRefiner:
    def __init__(self):
        pass

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
            if job_id:
                from app.api.scheduler import update_job_status
                await update_job_status(job_id, records_found=len(items))

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
                    if job_id:
                        await update_job_status(job_id, records_drafted=processed_count)
                    
                    logger.info(f"✅ Refined Knowledge Base Item: {item['id']}")
                    
                    # Rate limiting
                    await asyncio.sleep(0.5)
                    
                except Exception as e:
                    logger.error(f"Failed to refine item {item['id']}: {e}")
            
            return processed_count

        except Exception as e:
            logger.error(f"Knowledge Refiner Error: {e}")
            return 0

    async def analyze_abstract(self, abstract_text: str) -> Dict[str, Any]:
        """Gemini Flash를 사용하여 초록 분석 (Google SDK 직접 호출)"""
        
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
        
        full_prompt = f"{system_prompt}\n\nAbstract:\n{abstract_text}"
        
        try:
            # SDK 설정
            genai.configure(api_key=settings.GOOGLE_API_KEY)
            model = genai.GenerativeModel('gemini-2.0-flash-exp')
            
            logger.info(f"🚀 Requesting Gemini (Direct SDK) for PubMed abstract...")
            
            # 동기 호출을 비동기로 실행
            loop = asyncio.get_event_loop()
            response = await loop.run_in_executor(None, lambda: model.generate_content(full_prompt))
            
            content = response.text.strip()
            
            # 비용 추적
            usage = response.usage_metadata
            await cost_tracker.track_usage(
                "gemini-2.0-flash",
                usage.prompt_token_count,
                usage.candidates_token_count
            )
            
            # JSON 파싱 (json-repair 도입으로 백틱 공격 및 깨진 형식 방어)
            try:
                repaired_content = repair_json(content)
                return json.loads(repaired_content)
            except Exception as parse_e:
                logger.warning(f"Standard parsing failed for PubMed, trying manual strip: {parse_e}")
                if "```json" in content:
                    content = content.split("```json")[1].split("```")[0]
                elif "```" in content:
                    content = content.split("```")[1].split("```")[0]
                return json.loads(content.strip())
            
        except Exception as e:
            logger.error(f"AI Analysis Error: {e}")
            return {
                "summary": "Analysis Failed",
                "relevance_score": 0.0,
                "ai_reasoning": f"Error: {str(e)}"
            }

knowledge_refiner = KnowledgeBaseRefiner()
