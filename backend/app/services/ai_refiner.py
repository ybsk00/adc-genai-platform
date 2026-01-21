"""
AI Refiner Service
미정제 레코드를 LLM으로 분석하여 outcome_type, failure_reason 등을 추출
"""
import asyncio
import logging
import json
from typing import Dict, Any, Optional
from datetime import datetime

from app.core.supabase import supabase
from app.core.config import settings
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from app.services.cost_tracker import cost_tracker

logger = logging.getLogger(__name__)

class AIRefiner:
    def __init__(self):
        self.batch_size = 10  # 한 번에 처리할 레코드 수
        self.processed_count = 0
        self.error_count = 0

    async def _is_system_paused(self) -> bool:
        """시스템 일시정지 상태 확인"""
        try:
            res = supabase.table("system_config").select("value").eq("key", "AI_REFINER_STATUS").execute()
            return res.data[0]["value"] == "PAUSED" if res.data else False
        except:
            return False

    async def refine_single_record(self, record: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """단일 레코드 LLM 분석"""
        try:
            # 비용 한도 체크
            if await cost_tracker.is_over_limit():
                logger.warning("⚠️ Daily LLM cost limit reached. Skipping analysis.")
                return None
            
            llm = ChatOpenAI(
                model=settings.FAST_LLM,
                temperature=0,
                api_key=settings.OPENAI_API_KEY
            )
            
            # 원본 데이터에서 분석에 필요한 정보 추출
            properties = record.get("properties", {})
            raw_data = properties.get("raw_data", {})
            protocol = raw_data.get("protocolSection", {}) if raw_data else {}
            
            title = record.get("name", "No Title")
            description = properties.get("brief_summary", "")
            overall_status = properties.get("overall_status", "")
            why_stopped = properties.get("why_stopped", "")
            phase = properties.get("phase", "")
            
            # LLM 프롬프트
            system_prompt = """You are a Clinical Trial Analyst specializing in ADC (Antibody-Drug Conjugate) research.
Analyze the clinical trial data and extract structured information.

Output ONLY valid JSON in this exact format:
{
    "drug_name": "extracted drug/compound name or null",
    "target": "molecular target (e.g., HER2, TROP2) or null",
    "outcome_type": "Success|Failure|Ongoing|Unknown",
    "failure_reason": "reason if failed, null otherwise",
    "confidence": 0.0-1.0
}

Rules:
- outcome_type: "Success" if completed with positive results, "Failure" if terminated/withdrawn/negative, "Ongoing" if active, "Unknown" if unclear
- failure_reason: Only fill if outcome_type is "Failure"
- Be concise and accurate"""

            user_prompt = f"""Clinical Trial Analysis:
Title: {title}
Phase: {phase}
Status: {overall_status}
Why Stopped: {why_stopped}
Description: {description[:1000] if description else 'N/A'}"""

            prompt = ChatPromptTemplate.from_messages([
                ("system", system_prompt),
                ("user", user_prompt)
            ])
            
            chain = prompt | llm
            response = await chain.ainvoke({})
            content = response.content.strip()
            
            # 비용 추적
            usage = response.response_metadata.get('token_usage', {})
            await cost_tracker.track_usage(
                settings.FAST_LLM,
                usage.get('prompt_tokens', 0),
                usage.get('completion_tokens', 0)
            )
            
            # JSON 파싱
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0]
            elif "```" in content:
                content = content.split("```")[1].split("```")[0]
            
            analysis = json.loads(content)
            
            return {
                "drug_name": analysis.get("drug_name"),
                "target": analysis.get("target"),
                "outcome_type": analysis.get("outcome_type", "Unknown"),
                "failure_reason": analysis.get("failure_reason"),
                "ai_confidence": analysis.get("confidence", 0.5)
            }
        
        except Exception as e:
            logger.error(f"LLM Analysis Error for record {record.get('id')}: {e}")
            return None

    async def enrich_with_smiles(self, drug_name: Optional[str]) -> Optional[str]:
        """PubChem에서 SMILES 코드 조회"""
        if not drug_name:
            return None
        
        try:
            from app.services.chemical_resolver import chemical_resolver
            return chemical_resolver.fetch_verified_smiles(drug_name)
        except Exception as e:
            logger.error(f"SMILES lookup error for {drug_name}: {e}")
            return None

    async def process_pending_records(self, job_id: Optional[str] = None, max_records: int = 50):
        """
        미정제 레코드 배치 처리
        - ai_refined = false 인 레코드 조회
        - LLM 분석
        - SMILES 보강
        - DB 업데이트 (ai_refined = true)
        """
        from app.api.scheduler import update_job_status, is_cancelled
        
        # 시스템 일시정지 체크
        if await self._is_system_paused():
            logger.info("⏸️ System is PAUSED. Skipping AI Refiner.")
            if job_id:
                await update_job_status(job_id, status="skipped", errors=["System is paused"])
            return
        
        if job_id:
            await update_job_status(job_id, status="running")
        
        logger.info("🧹 [AI Refiner] Checking for pending records...")
        
        try:
            # 미정제 레코드 조회
            response = supabase.table("golden_set_library")\
                .select("*")\
                .eq("ai_refined", False)\
                .limit(max_records)\
                .execute()
            
            pending_items = response.data
            
            if not pending_items:
                logger.info("✨ No pending items. Golden Set is clean.")
                if job_id:
                    await update_job_status(
                        job_id, 
                        status="completed", 
                        records_found=0,
                        completed_at=datetime.utcnow().isoformat()
                    )
                return
            
            if job_id:
                await update_job_status(job_id, records_found=len(pending_items))
            
            logger.info(f"🔍 [AI Refiner] Processing {len(pending_items)} items...")
            
            refined_count = 0
            error_count = 0
            
            for item in pending_items:
                # 중단 요청 체크
                if job_id and await is_cancelled(job_id):
                    logger.info("Refiner cancelled by user")
                    await update_job_status(job_id, status="stopped")
                    return
                
                try:
                    # LLM 분석
                    analysis = await self.refine_single_record(item)
                    
                    if analysis:
                        # SMILES 보강
                        smiles = await self.enrich_with_smiles(analysis.get("drug_name"))
                        
                        # 기존 properties에 AI 분석 결과 추가
                        updated_properties = item.get("properties", {})
                        updated_properties["ai_analysis"] = analysis
                        # raw_data는 너무 크므로 제거
                        if "raw_data" in updated_properties:
                            del updated_properties["raw_data"]
                        
                        # DB 업데이트
                        update_payload = {
                            "name": analysis.get("drug_name") or item.get("name"),
                            "outcome_type": analysis.get("outcome_type", "Unknown"),
                            "failure_reason": analysis.get("failure_reason"),
                            "smiles_code": smiles,
                            "ai_refined": True,
                            "properties": updated_properties
                        }
                        
                        supabase.table("golden_set_library")\
                            .update(update_payload)\
                            .eq("id", item["id"])\
                            .execute()
                        
                        refined_count += 1
                        logger.info(f"✅ Refined: {item.get('name', 'Unknown')[:50]}...")
                    else:
                        # 분석 실패 시 에러 기록
                        supabase.table("golden_set_library")\
                            .update({
                                "processing_error": "LLM analysis failed",
                                "ai_refined": True  # 재시도 방지를 위해 true로 설정
                            })\
                            .eq("id", item["id"])\
                            .execute()
                        error_count += 1
                    
                    # 진행률 업데이트
                    if job_id:
                        await update_job_status(job_id, records_drafted=refined_count)
                    
                    # Rate limiting
                    await asyncio.sleep(0.5)
                
                except Exception as e:
                    logger.error(f"Record processing error: {e}")
                    error_count += 1
            
            logger.info(f"🎉 Refiner Complete! Refined: {refined_count}, Errors: {error_count}")
            
            if job_id:
                await update_job_status(
                    job_id,
                    status="completed",
                    records_drafted=refined_count,
                    completed_at=datetime.utcnow().isoformat(),
                    errors=[f"Errors: {error_count}"] if error_count > 0 else []
                )
        
        except Exception as e:
            logger.error(f"AI Refiner Error: {e}")
            if job_id:
                await update_job_status(job_id, status="failed", errors=[str(e)])

# 싱글톤 인스턴스
ai_refiner = AIRefiner()
