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
logger.setLevel(logging.DEBUG)

class AIRefiner:
    # ADC 관련 키워드 (Pre-filtering용)
    ADC_KEYWORDS = [
        "antibody-drug conjugate", "adc", "immunoconjugate",
        "trastuzumab", "vedotin", "emtansine", "ozogamicin", 
        "deruxtecan", "govitecan", "mertansine", "ravtansine",
        "duocarmycin", "maytansine", "auristatin", "calicheamicin",
        "her2", "trop2", "cd19", "cd22", "cd33", "cd79", "bcma",
        "nectin-4", "folate receptor", "egfr", "psma"
    ]
    
    # Non-cancer 질환 제외 키워드
    EXCLUDE_KEYWORDS = [
        "alzheimer", "diabetes", "parkinson", "arthritis",
        "hypertension", "cardiovascular", "obesity", "asthma",
        "copd", "depression", "anxiety", "schizophrenia"
    ]
    
    # Gemini Safety Settings (의학 용어 차단 해제)
    SAFETY_SETTINGS = [
        {"category": "HARM_CATEGORY_HATE_SPEECH", "threshold": "BLOCK_NONE"},
        {"category": "HARM_CATEGORY_HARASSMENT", "threshold": "BLOCK_NONE"},
        {"category": "HARM_CATEGORY_SEXUALLY_EXPLICIT", "threshold": "BLOCK_NONE"},
        {"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}
    ]
    
    def __init__(self):
        logger.info("🔥 [AI Refiner] Service Initialized (Version: 2026-01-24-1750)")
        self.batch_size = 10  # 한 번에 처리할 레코드 수
        self.processed_count = 0
        self.error_count = 0
        self.semaphore = asyncio.Semaphore(10) # 동시성 제어 (최대 10개)
    
    def is_adc_relevant(self, record: Dict[str, Any]) -> bool:
        """ADC 관련 데이터인지 Pre-filter 체크"""
        properties = record.get("properties", {})
        title = (record.get("name", "") or record.get("product_name", "") or "").lower()
        description = (properties.get("brief_summary", "") or record.get("summary", "") or "").lower()
        indication = (properties.get("indication", "") or "").lower()
        moa = (properties.get("mechanism_of_action", "") or "").lower()
        
        combined_text = f"{title} {description} {indication} {moa}"
        
        # 1. Non-cancer 질환 제외
        for exclude in self.EXCLUDE_KEYWORDS:
            if exclude in combined_text:
                logger.info(f"⛔ Pre-filter SKIP (Non-cancer): {record.get('name', 'Unknown')[:40]} - contains '{exclude}'")
                return False
        
        # 2. ADC 키워드 필수 포함 체크
        for keyword in self.ADC_KEYWORDS:
            if keyword in combined_text:
                return True
        
        logger.info(f"⛔ Pre-filter SKIP (No ADC keyword): {record.get('name', 'Unknown')[:40]}")
        return False

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
            logger.info(f"⚡ [AI Refiner] Starting analysis for record {record.get('id')} ({record.get('source_name', 'Unknown')})")

            # 비용 한도 체크
            if await cost_tracker.is_over_limit():
                logger.warning("⚠️ Daily LLM cost limit reached. Skipping analysis.")
                return None
            
            # SDK 설정
            genai.configure(api_key=settings.GOOGLE_API_KEY)
            model_id = settings.GEMINI_MODEL_ID or 'gemini-2.0-flash'
            model = genai.GenerativeModel(model_id)
            
            # 원본 데이터에서 분석에 필요한 정보 추출
            properties = record.get("properties", {})
            title = record.get("name", "No Title")
            
            # enrichment_source에 따라 텍스트 추출 분기
            source = record.get("enrichment_source")
            
            if source == "open_fda_api":
                # OpenFDA: fda_label에서 데이터 추출 (다중 경로 시도)
                fda_label = properties.get("fda_label", {})
                
                # 1. Indication 추출 (여러 필드 순차 시도)
                description = (
                    properties.get("indication", "") or
                    fda_label.get("indication", "") or 
                    fda_label.get("indications_and_usage", "") or 
                    properties.get("indications_and_usage", "") or
                    fda_label.get("description", "") or 
                    properties.get("description", "")
                )
                
                # 2. Mechanism of Action 추출 (타겟 정보 포함)
                moa = (
                    properties.get("mechanism_of_action", "") or
                    fda_label.get("mechanism_of_action", "")
                )
                
                # 3. Boxed Warning 추출
                boxed_warning = (
                    properties.get("boxed_warning", "") or
                    fda_label.get("boxed_warning", "") or 
                    fda_label.get("warnings", "")
                )
                
                # 4. Generic Name 추출
                generic_name = (
                    properties.get("generic_name", "") or
                    fda_label.get("generic_name", "")
                )
                
                # 🔍 디버그 로그: 추출된 텍스트 확인
                logger.info(f"📋 [OpenFDA Parse] Drug: {title}")
                logger.info(f"   - Description length: {len(description)} chars")
                logger.info(f"   - MoA length: {len(moa)} chars")
                logger.info(f"   - Generic Name: {generic_name or 'N/A'}")
                
                # 빈 텍스트 경고
                if not description and not moa:
                    logger.warning(f"⚠️ No text found for {title}! Properties keys: {list(properties.keys())}")
                    if fda_label:
                        logger.warning(f"   fda_label keys: {list(fda_label.keys())}")
            elif record.get("source_name") in ["Ambeed", "Creative Biolabs"]:
                # 상용 시약 데이터
                description = record.get("summary", "") or record.get("product_name", "")
                moa = ""
                boxed_warning = ""
                generic_name = record.get("cas_number", "")
            else:
                # Clinical Trials: brief_summary 사용
                description = properties.get("brief_summary", "")
                moa = ""
                boxed_warning = ""
                generic_name = ""
            
            overall_status = properties.get("overall_status", "")
            why_stopped = properties.get("why_stopped", "")
            phase = properties.get("phase", "")
            
            if source == "open_fda_api":
                # OpenFDA 전용 프롬프트 (Indication 우선 타겟 추출 + ADC 키워드 감지)
                system_prompt = """You are a Pharmaceutical Regulatory Affairs Specialist analyzing FDA Drug Labels.

CRITICAL INSTRUCTIONS:
1. Target Extraction Priority:
   - FIRST check the Indication text for molecular targets (e.g., HER2, CD19, FRα, FR-alpha, TROP2, CD22, CD33, Nectin-4)
   - If not found in Indication, check Mechanism of Action
   - Common patterns: "[drug] targets [protein]", "[protein]-positive", "[protein]-directed"

2. ADC Detection for Relevance Score:
   - If the drug name ends with "-mab" or "-tansine" or contains "conjugate", set relevance_score = 1.0
   - If text mentions "antibody-drug conjugate", "immunoconjugate", "ADC", set relevance_score = 1.0
   - Otherwise, set relevance_score based on how related to ADC research (0.0-1.0)

Output ONLY valid JSON:
{
    "drug_name": "extracted drug name",
    "target": "molecular target (e.g., HER2, CD19, FRα, TROP2) - MUST extract from Indication if available",
    "outcome_type": "Success",
    "approval_status": "Approved",
    "boxed_warning": "Summary of Boxed Warning or 'None'",
    "indication": "Primary cancer/disease indication",
    "binding_affinity": "extracted Kd value (e.g. 1.2 nM) or null",
    "isotype": "extracted isotype (e.g. IgG1) or null",
    "host_species": "extracted host (e.g. Human) or null",
    "orr_pct": "ORR percentage value (number only) or null",
    "os_months": "OS in months (number only) or null",
    "pfs_months": "PFS in months (number only) or null",
    "relevance_score": 0.0-1.0,
    "confidence": 0.0-1.0
}
"""
                full_prompt = f"""{system_prompt}

FDA Label Data:
Name: {title}
Generic Name: {generic_name}
Indication: {description[:800] if description else "N/A"}
Mechanism of Action: {moa[:500] if moa else "N/A"}
Boxed Warning: {boxed_warning[:200] if boxed_warning else "N/A"}
"""
                # 🔍 디버그: 최종 프롬프트 길이 로그
                logger.info(f"📤 Gemini Prompt Length: {len(full_prompt)} chars")
            elif record.get("source_name") in ["Ambeed", "Creative Biolabs"]:
                # 상용 시약 전용 프롬프트
                system_prompt = """You are a Bio-Chemical Specialist analyzing ADC Reagents.
Analyze the product data and classify it.

Output ONLY valid JSON:
{
    "drug_name": "product name",
    "target": "molecular target (e.g., HER2, TROP2) if applicable, else null",
    "category": "Target|Payload|Linker|Drug-Linker|Conjugate|Other",
    "outcome_type": "Success",
    "relevance_score": 0.0-1.0,
    "confidence": 0.0-1.0
}
"""
                full_prompt = f"""{system_prompt}

Product Data:
Name: {record.get('product_name')}
Category: {record.get('category')}
CAS: {record.get('cas_number')}
Summary: {description}
"""
            else:
                # 기존 Clinical Trials 프롬프트
                system_prompt = """You are a Clinical Trial Analyst specializing in ADC (Antibody-Drug Conjugate) research.
Analyze the clinical trial data and extract structured information.

Output ONLY valid JSON in this exact format:
{
    "drug_name": "extracted drug/compound name or null",
    "target": "molecular target (e.g., HER2, TROP2) or null",
    "outcome_type": "Success|Failure|Ongoing|Unknown",
    "failure_reason": "reason if failed, null otherwise",
    "binding_affinity": "extracted Kd value or null",
    "isotype": "extracted isotype or null",
    "host_species": "extracted host or null",
    "orr_pct": "ORR percentage or null",
    "os_months": "OS in months or null",
    "pfs_months": "PFS in months or null",
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

            logger.info(f"🚀 Requesting Gemini (Direct SDK) for record {record.get('id')} ({source})...")
            
            # 동기 호출을 비동기로 실행 (Safety Settings 적용)
            loop = asyncio.get_event_loop()
            response = await loop.run_in_executor(None, lambda: model.generate_content(
                full_prompt,
                safety_settings=self.SAFETY_SETTINGS
            ))
            
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
            
            return {
                "drug_name": analysis.get("drug_name"),
                "target": analysis.get("target"),
                "outcome_type": analysis.get("outcome_type", "Unknown"),
                "failure_reason": analysis.get("failure_reason"),
                "binding_affinity": analysis.get("binding_affinity"),
                "isotype": analysis.get("isotype"),
                "host_species": analysis.get("host_species"),
                "orr_pct": analysis.get("orr_pct"),
                "os_months": analysis.get("os_months"),
                "pfs_months": analysis.get("pfs_months"),
                "ai_confidence": analysis.get("confidence", 0.5),
                "relevance_score": analysis.get("relevance_score", 0.0),
                "boxed_warning": analysis.get("boxed_warning"), # OpenFDA specific
                "indication": analysis.get("indication"), # OpenFDA specific
                "generic_name": generic_name # Pass generic name back
            }
        
        except Exception as e:
            logger.error(f"LLM Analysis Error for record {record.get('id')}: {e}")
            # 에러 로그를 sync_jobs에 남기기 (비동기로 시도)
            try:
                error_data = {
                    "status": "failed",
                    "errors": [f"LLM Error for {record.get('id')}: {str(e)}"],
                    "completed_at": datetime.utcnow().isoformat()
                }
                # job_id가 없으면 새로 생성하거나 로그용 job_id 사용
                supabase.table("sync_jobs").insert({
                    "status": "failed",
                    "errors": [f"LLM Error: {str(e)}"],
                    "records_found": 1,
                    "records_drafted": 0
                }).execute()
            except:
                pass
            return {"error": str(e)}

    async def enrich_with_pubchem(self, drug_name: Optional[str], generic_name: Optional[str] = None) -> Dict[str, Any]:
        """2단계: PubChem 화학 구조 자동 매핑 (Fuzzy Matching 포함)"""
        if not drug_name and not generic_name:
            return {}
        
        try:
            import pubchempy as pcp
            from app.services.chemical_resolver import chemical_resolver
            
            loop = asyncio.get_event_loop()
            
            def fetch_pubchem():
                # Chemical Resolver의 향상된 fetch_safe_smiles 사용 (Generic Fallback 포함)
                result = chemical_resolver.fetch_safe_smiles(drug_name, generic_name)
                
                if result["smiles"]:
                    return {
                        "smiles_code": result["smiles"],
                        "canonical_smiles": result["smiles"], # fetch_safe_smiles returns isomeric/canonical
                        "molecular_weight": result["mw"],
                        "enrichment_source": "PubChem"
                    }
                return None

            result = await loop.run_in_executor(None, fetch_pubchem)
            return result
            
        except Exception as e:
            logger.error(f"PubChem lookup error for {drug_name}: {e}")
            return {"error": "PubChem Error"}
            
        except Exception as e:
            logger.error(f"PubChem lookup error for {drug_name}: {e}")
            return {"error": "PubChem Error"}

    async def generate_smiles_with_ai(self, drug_name: str) -> Dict[str, Any]:
        """Gemini Fallback: 화학자 페르소나로 SMILES 생성"""
        try:
            genai.configure(api_key=settings.GOOGLE_API_KEY)
            model = genai.GenerativeModel('gemini-2.5-flash')
            
            prompt = f"""You are an expert computational chemist. 
Generate the Canonical SMILES for the drug "{drug_name}".
If the exact structure is unknown, infer it from the most common derivative or similar structure.
Output ONLY the SMILES string. Do not include any explanation or markdown."""

            loop = asyncio.get_event_loop()
            response = await loop.run_in_executor(None, lambda: model.generate_content(prompt))
            smiles = response.text.strip().replace("```", "").strip()
            
            # Sanity Check with RDKit
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return {
                    "smiles_code": smiles,
                    "canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True), # 정규화
                    "enrichment_source": "AI-Generated"
                }
            else:
                logger.warning(f"⚠️ AI generated invalid SMILES for {drug_name}: {smiles}")
                return {"error": "Invalid SMILES generated"}
                
        except Exception as e:
            logger.error(f"AI SMILES generation failed: {e}")
            return {"error": str(e)}

    async def process_pending_records(self, job_id: Optional[str] = None, max_records: int = 50, source_filter: Optional[str] = None):
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
        
        logger.info(f"🧹 [AI Refiner] Checking for pending records (Source: {source_filter or 'All'})...")
        
        try:
            # 미정제 레코드 조회
            query = supabase.table("golden_set_library")\
                .select("*")\
                .eq("ai_refined", False)
            
            if source_filter:
                query = query.eq("enrichment_source", source_filter)
                
            response = query.limit(max_records).execute()
            
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
                
                # ⚡ Pre-filter: ADC 관련 데이터만 처리 (비용 절감)
                if not self.is_adc_relevant(item):
                    # 노이즈 데이터는 Skip 처리하고 마킹
                    supabase.table("golden_set_library")\
                        .update({
                            "ai_refined": True,
                            "rag_status": "excluded",
                            "processing_error": "Pre-filter: Not ADC related",
                            "relevance_score": 0.0
                        })\
                        .eq("id", item["id"])\
                        .execute()
                    continue
                
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
                    
                    async with self.semaphore: # 세마포어 적용
                        if analysis and "error" not in analysis:
                            drug_name = analysis.get("drug_name") or existing_drug_name
                            relevance_score = analysis.get("relevance_score", 0.0)
                            generic_name = analysis.get("generic_name") or item.get("properties", {}).get("generic_name")
                            
                            # 2️⃣ 화학 구조 매핑 (관련성 점수와 무관하게 이름이 있으면 시도)
                            # 사용자 요청: "브랜드명으로 SMILES 못 찾으면 성분명으로 끝까지 찾아내는 로직이 핵심"
                            pubchem_data = None
                            processing_error = None
                            
                            # SMILES 조회 조건 완화: 이름만 있으면 무조건 시도
                            if drug_name or generic_name:
                                if existing_smiles:
                                    logger.info(f"⏩ PubChem Skip: {drug_name[:30]} (SMILES exists)")
                                    pubchem_data = {"smiles_code": existing_smiles, "enrichment_source": "Existing"}
                                else:
                                    # 2-1. PubChem Lookup (with Generic Fallback)
                                    logger.info(f"🔬 PubChem lookup: {drug_name} (Generic: {generic_name})")
                                    pubchem_data = await self.enrich_with_pubchem(drug_name, generic_name)
                                    
                                    # 2-2. Fallback to AI (Only if PubChem failed completely)
                                    if not pubchem_data or "error" in pubchem_data:
                                        target_name = drug_name or generic_name
                                        logger.info(f"🧪 AI Fallback: Generating SMILES for {target_name}")
                                        pubchem_data = await self.generate_smiles_with_ai(target_name)
                                        
                                        if not pubchem_data or "error" in pubchem_data:
                                            processing_error = "SMILES Not Found (PubChem & AI Failed)"
                                            logger.warning(f"⚠️ All methods failed for: {target_name}")
                            else:
                                logger.warning(f"⚠️ No drug name found for SMILES lookup: {item.get('id')}")
                            
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
                                "rag_status": "processed", 
                                "processing_error": processing_error,
                                "properties": updated_properties
                            }
                            
                            # PubChem/AI 데이터 병합
                            if pubchem_data and "error" not in pubchem_data:
                                update_payload.update(pubchem_data)
                            
                            supabase.table("golden_set_library")\
                                .update(update_payload)\
                                .eq("id", item["id"])\
                                .execute()
                            
                            refined_count += 1
                            source = pubchem_data.get("enrichment_source", "None") if pubchem_data else "None"
                            logger.info(f"✅ Refined: {drug_name[:30]} (Score: {relevance_score}, Source: {source})")
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
