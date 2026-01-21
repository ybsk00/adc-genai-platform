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
import google.generativeai as genai
from app.services.cost_tracker import cost_tracker
from json_repair import repair_json

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
        """단일 레코드 LLM 분석 (Google SDK 직접 호출)"""
        try:
            # 비용 한도 체크
            if await cost_tracker.is_over_limit():
                logger.warning("⚠️ Daily LLM cost limit reached. Skipping analysis.")
                return None
            
            # SDK 설정
            genai.configure(api_key=settings.GOOGLE_API_KEY)
            model = genai.GenerativeModel('gemini-2.0-flash-exp')
            
            # 원본 데이터에서 분석에 필요한 정보 추출
            properties = record.get("properties", {})
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
    "relevance_score": 0.0-1.0 (relevance to ADC research),
    "confidence": 0.0-1.0
}

Rules:
- outcome_type: "Success" if completed with positive results, "Failure" if terminated/withdrawn/negative, "Ongoing" if active, "Unknown" if unclear
- failure_reason: Only fill if outcome_type is "Failure"
- Be concise and accurate

IMPORTANT: Return ONLY raw JSON. Do not use markdown formatting like ```json ... ```.
"""

            full_prompt = f"""{system_prompt}

Clinical Trial Analysis:
Title: {title}
Phase: {phase}
Status: {overall_status}
Why Stopped: {why_stopped}
Description: {description[:1000] if description else 'N/A'}"""

            logger.info(f"🚀 Requesting Gemini (Direct SDK) for record {record.get('id')}...")
            
            # 동기 호출을 비동기로 실행
            loop = asyncio.get_event_loop()
            response = await loop.run_in_executor(None, lambda: model.generate_content(full_prompt))
            
            content = response.text.strip()
            
            # 비용 추적 (Gemini 2.0 Flash 대략적 토큰 계산 - SDK에서 직접 가져오기 어려울 경우 대비)
            # 실제로는 response.usage_metadata에 있음
            usage = response.usage_metadata
            await cost_tracker.track_usage(
                "gemini-2.0-flash",
                usage.prompt_token_count,
                usage.candidates_token_count
            )
            
            # JSON 파싱 (json-repair 도입으로 백틱 공격 및 깨진 형식 방어)
            try:
                repaired_content = repair_json(content)
                analysis = json.loads(repaired_content)
            except Exception as parse_e:
                logger.warning(f"Standard parsing failed, trying manual strip: {parse_e}")
                if "```json" in content:
                    content = content.split("```json")[1].split("```")[0]
                elif "```" in content:
                    content = content.split("```")[1].split("```")[0]
                analysis = json.loads(content.strip())
            
            # 3단계 파이프라인: 1. 스마트 필터링 & 추출
            return {
                "drug_name": analysis.get("drug_name"),
                "target": analysis.get("target"),
                "outcome_type": analysis.get("outcome_type", "Unknown"),
                "failure_reason": analysis.get("failure_reason"),
                "ai_confidence": analysis.get("confidence", 0.5),
                "relevance_score": analysis.get("relevance_score", 0.0) # 프롬프트에 추가 필요
            }
        
        except Exception as e:
            logger.error(f"LLM Analysis Error for record {record.get('id')}: {e}")
            return {"error": str(e)}

    async def enrich_with_pubchem(self, drug_name: Optional[str]) -> Dict[str, Any]:
        """2단계: PubChem 화학 구조 자동 매핑"""
        if not drug_name:
            return {}
        
        try:
            # pubchempy는 동기 라이브러리이므로 비동기로 실행
            import pubchempy as pcp
            loop = asyncio.get_event_loop()
            
            def fetch_pubchem():
                compounds = pcp.get_compounds(drug_name, 'name')
                if compounds:
                    c = compounds[0]
                    return {
                        "smiles_code": c.isomeric_smiles,
                        "canonical_smiles": c.canonical_smiles,
                        "molecular_weight": float(c.molecular_weight) if c.molecular_weight else None
                    }
                return None

            result = await loop.run_in_executor(None, fetch_pubchem)
            return result
            
        except Exception as e:
            logger.error(f"PubChem lookup error for {drug_name}: {e}")
            return {"error": "PubChem Error"}

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
                    # ============ Smart Skip Logic ============
                    existing_drug_name = item.get("name")
                    existing_outcome = item.get("outcome_type")
                    existing_smiles = item.get("smiles_code")
                    
                    # 1️⃣ LLM 분석 (drug_name 또는 outcome_type이 없을 때만)
                    if not existing_drug_name or existing_drug_name == "Unknown" or not existing_outcome:
                        logger.info(f"🤖 LLM analyzing: {item.get('id', 'unknown')[:20]}...")
                        analysis = await self.refine_single_record(item)
                    else:
                        logger.info(f"⏩ LLM Skip: {existing_drug_name[:30]} (already extracted)")
                        analysis = {
                            "drug_name": existing_drug_name,
                            "outcome_type": existing_outcome,
                            "failure_reason": item.get("failure_reason"),
                            "target": item.get("properties", {}).get("target")
                        }
                    
                    if analysis and "error" not in analysis:
                        drug_name = analysis.get("drug_name") or existing_drug_name
                        relevance_score = analysis.get("relevance_score", 0.0)
                        
                        # 2️⃣ 화학 구조 매핑 (관련성 높을 때만)
                        pubchem_data = None
                        processing_error = None
                        
                        if relevance_score > 0.5 and drug_name:
                            if existing_smiles:
                                logger.info(f"⏩ PubChem Skip: {drug_name[:30]} (SMILES exists)")
                                pubchem_data = {"smiles_code": existing_smiles}
                            else:
                                logger.info(f"🔬 PubChem lookup: {drug_name}")
                                pubchem_data = await self.enrich_with_pubchem(drug_name)
                                
                                if not pubchem_data:
                                    # [실패 대응] PubChem 조회 실패 시 흔적 남기기
                                    processing_error = "PubChem Not Found"
                                    logger.warning(f"⚠️ PubChem Not Found for: {drug_name}")
                        
                        # 기존 properties에 AI 분석 결과 추가
                        updated_properties = item.get("properties", {})
                        updated_properties["ai_analysis"] = analysis
                        if "raw_data" in updated_properties:
                            del updated_properties["raw_data"]
                        
                        # 3️⃣ DB 상태 업데이트
                        update_payload = {
                            "name": drug_name or item.get("name"),
                            "outcome_type": analysis.get("outcome_type", "Unknown"),
                            "failure_reason": analysis.get("failure_reason"),
                            "relevance_score": relevance_score,
                            "ai_refined": True,
                            "rag_status": "processed", # 성공 처리
                            "processing_error": processing_error, # PubChem 실패 시 기록
                            "properties": updated_properties
                        }
                        
                        # PubChem 데이터 병합
                        if pubchem_data and "error" not in pubchem_data:
                            update_payload.update(pubchem_data)
                        
                        supabase.table("golden_set_library")\
                            .update(update_payload)\
                            .eq("id", item["id"])\
                            .execute()
                        
                        refined_count += 1
                        logger.info(f"✅ Refined: {drug_name[:50] if drug_name else 'Unknown'} (Score: {relevance_score})")
                    else:
                        # 분석 실패 시 상세 에러 기록
                        error_msg = analysis.get("error") if analysis else "Unknown LLM Error"
                        supabase.table("golden_set_library")\
                            .update({
                                "processing_error": error_msg,
                                "ai_refined": True,
                                "rag_status": "failed"
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
                    logger.error(f"Record processing error for {item.get('id')}: {e}")
                    # 에러 발생 시에도 처리된 것으로 표시하여 무한 루프 방지
                    try:
                        supabase.table("golden_set_library")\
                            .update({
                                "processing_error": f"Processing failed: {str(e)}",
                                "ai_refined": True, # 재시도 방지
                                "outcome_type": "Failure",
                                "failure_reason": "System Error"
                            })\
                            .eq("id", item["id"])\
                            .execute()
                    except Exception as db_e:
                        logger.error(f"Failed to update error status: {db_e}")
                    
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
