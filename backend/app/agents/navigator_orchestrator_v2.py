"""
One-Click ADC Navigator Orchestrator V2
에이전트 실제 호출 + DesignSessionState 통합 버전

FIXED:
- Librarian, Alchemist, Coder, Healer, Auditor 실제 execute() 호출
- DesignSessionState 공유 메모리 사용
- 실시간 연산 로그 스트리밍
- Healer 자동 재시도 연동
- 데이터 출처 추적 (PMID, 에이전트별 기록)
- f-string 인젝션 방지
"""

import logging
import json
import uuid
import asyncio
import re
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum

from app.core.supabase import get_supabase_client
from app.core.websocket_hub import websocket_hub
from app.core.config import settings
from app.core.gemini import get_gemini_model

# FIXED: Import actual agents
from app.agents.librarian import LibrarianAgent
from app.agents.alchemist import AlchemistAgent
from app.agents.coder import CoderAgent
from app.agents.healer import HealerAgent
from app.agents.auditor import AuditorAgent

# FIXED: Import shared state
from app.agents.design_state import (
    DesignSessionState, 
    create_initial_state,
    CandidateStructure,
    CalculatedMetrics,
    ValidationFlags
)
from app.agents.base_agent import AgentOutput

logger = logging.getLogger(__name__)


# ============================================================================
# Navigator Result Types
# ============================================================================

@dataclass
class NavigatorAgentLog:
    """에이전트 실행 로그 (실시간 스트리밍용)"""
    timestamp: str
    agent_name: str
    step: int
    status: str  # started, running, completed, error, healing
    message: str
    reasoning: Optional[str] = None
    data_source: Optional[str] = None  # 데이터 출처 추적
    pmid_references: List[str] = field(default_factory=list)
    confidence_score: Optional[float] = None


@dataclass
class NavigatorResultV2:
    """Navigator V2 결과"""
    session_id: str
    disease_name: str
    antibody_candidates: List[Dict[str, Any]]
    golden_combination: Dict[str, Any]
    calculated_metrics: Dict[str, Any]
    physics_verified: bool
    virtual_trial: Dict[str, Any]
    agent_logs: List[NavigatorAgentLog]  # 실시간 연산 로그
    data_lineage: Dict[str, Any]  # 데이터 출처 추적
    execution_time_seconds: float = 0.0
    warnings: List[str] = field(default_factory=list)
    critical_errors: List[str] = field(default_factory=list)  # Fail-Fast: 사용자에게 노출할 치명적 오류


# ============================================================================
# Navigator Orchestrator V2 - Agent Integration
# ============================================================================

class NavigatorOrchestratorV2:
    """
    One-Click ADC Navigator V2
    
    FIXED: 실제 에이전트 호출 + DesignSessionState 공유
    """

    def __init__(self):
        self.supabase = get_supabase_client()
        
        # FIXED: 실제 에이전트 인스턴스 생성
        self.librarian = LibrarianAgent()
        self.alchemist = AlchemistAgent()
        self.coder = CoderAgent()
        self.healer = HealerAgent()
        self.auditor = AuditorAgent()
        
        self.agent_logs: List[NavigatorAgentLog] = []

    async def run_navigator_v2(
        self,
        disease_name: str,
        session_id: Optional[str] = None,
        user_id: Optional[str] = None,
        target_protein: Optional[str] = None,
        selected_antibody_id: Optional[str] = None
    ) -> NavigatorResultV2:
        """
        Navigator V2 실행 - 실제 에이전트 호출 + Selection Gate
        """
        if not session_id:
            session_id = str(uuid.uuid4())

        start_time = datetime.utcnow()
        warnings = []

        # FIXED: 세션별 agent_logs 초기화 (싱글톤 공유 문제 해결)
        self.agent_logs = []
        
        # FIXED: DesignSessionState 생성 (사용자 선택 타겟 반영)
        effective_target = target_protein or disease_name
        state = create_initial_state(
            session_id=session_id,
            user_id=user_id or "anonymous",
            session_type="navigator",
            tier="premium",
            target_antigen=effective_target,
            target_indication=disease_name,
            requested_dar=4,
            linker_preference="cleavable",
            design_goal=f"Design optimal ADC targeting {effective_target} for {disease_name}"
        )
        
        # 초기 상태 설정
        state["status"] = "running"
        state["step"] = 0
        
        target_msg = f" (target: {target_protein})" if target_protein else ""
        await self._log_agent_event(session_id, "orchestrator", 0, "started",
                                    f"Starting Navigator V2 for {disease_name}{target_msg}")

        try:
            # ═══════════════════════════════════════════════════════════
            # Step 1: Librarian - Target & Antibody Match
            # ═══════════════════════════════════════════════════════════
            await self._log_agent_event(session_id, "librarian", 1, "started",
                                        "Searching knowledge base for disease information...")
            
            await self._broadcast_step(session_id, 1, "Librarian: Disease analysis...")
            
            # FIXED: Librarian 실제 호출
            librarian_result = await self.librarian.execute(state)
            
            if not librarian_result.success:
                raise ValueError(f"Librarian failed: {librarian_result.error}")
            
            # 결과 처리
            antibody_candidates = self._extract_antibodies_from_librarian(librarian_result, disease_name)
            pmid_refs = librarian_result.referenced_pmids or []
            
            await self._log_agent_event(
                session_id, "librarian", 1, "completed",
                f"Found {len(antibody_candidates)} antibody candidates",
                reasoning=librarian_result.reasoning,
                data_source=librarian_result.data.get("normalized_target"),
                pmid_references=pmid_refs,
                confidence_score=librarian_result.confidence_score
            )

            # [gpNMB / HIGH-RISK TARGET WARNING]
            for ab in antibody_candidates:
                target = ab.get("target_protein", "")
                if target in self.HIGH_RISK_TARGETS:
                    warning_msg = self.HIGH_RISK_TARGETS[target]
                    warnings.append(f"HIGH-RISK TARGET ({target}): {warning_msg}")
                    await self._log_agent_event(
                        session_id, "orchestrator", 1, "warning",
                        f"Target {target} has failed clinical trials: {warning_msg}",
                        confidence_score=0.3
                    )

            # 사용자 선택 타겟이 있으면 해당 타겟 항체를 우선 배치
            if target_protein and antibody_candidates:
                # 선택된 타겟에 매칭되는 항체를 맨 앞으로
                matched = [ab for ab in antibody_candidates if ab.get("target_protein") == target_protein]
                others = [ab for ab in antibody_candidates if ab.get("target_protein") != target_protein]
                antibody_candidates = matched + others

            # Auto-select top antibody
            if antibody_candidates:
                top_ab = antibody_candidates[0]
                state["target_antigen"] = top_ab.get("target_protein", target_protein or disease_name)
                await self._log_agent_event(
                    session_id, "orchestrator", 1, "auto_select",
                    f"Selected antibody: {top_ab.get('name')} (target: {top_ab.get('target_protein')})"
                )

                # 항체 후보 저장
                try:
                    self.supabase.table("navigator_sessions").update({
                        "antibody_candidates": antibody_candidates,
                        "primary_target": top_ab.get("target_protein"),
                        "current_step": 1
                    }).eq("id", session_id).execute()
                except Exception as e:
                    logger.warning(f"[navigator-v2] Failed to save antibody candidates: {e}")
            
            # ═══════════════════════════════════════════════════════════
            # Step 2: Alchemist - Golden Combination
            # ═══════════════════════════════════════════════════════════
            await self._log_agent_event(session_id, "alchemist", 2, "started",
                                        "Designing ADC structure with Golden Set references...")
            
            await self._broadcast_step(session_id, 2, "Alchemist: Structure design...")
            
            # FIXED: Alchemist 실제 호출
            # Librarian 결과를 state에 반영
            state["pmid_references"] = pmid_refs
            state["golden_set_references"] = librarian_result.data.get("golden_set_references", [])
            
            alchemist_result = await self.alchemist.execute(state)
            
            if not alchemist_result.success:
                raise ValueError(f"Alchemist failed: {alchemist_result.error}")
            
            # 결과 반영
            candidates = alchemist_result.data.get("candidates", [])
            if candidates:
                state["candidates"] = candidates
                state["current_smiles"] = candidates[0].get("smiles", "")
            
            await self._log_agent_event(
                session_id, "alchemist", 2, "completed",
                f"Generated {len(candidates)} candidate structures",
                reasoning=alchemist_result.reasoning,
                data_source=f"golden_set:{len(alchemist_result.data.get('golden_set_refs', []))}, "
                           f"reagents:{len(alchemist_result.data.get('reagent_refs', []))}",
                confidence_score=alchemist_result.confidence_score
            )

            # [PBD / HIGH-TOXICITY PAYLOAD GUARDRAIL]
            payload_class = alchemist_result.data.get("payload_class", "")
            if not payload_class and candidates:
                payload_class = candidates[0].get("payload_class", "")
            for toxic_payload, risk_msg in self.HIGH_TOXICITY_PAYLOADS.items():
                if toxic_payload.lower() in payload_class.lower():
                    warnings.append(f"HIGH-TOXICITY PAYLOAD ({toxic_payload}): {risk_msg}")
                    await self._log_agent_event(
                        session_id, "orchestrator", 2, "warning",
                        f"Payload class {toxic_payload} flagged: {risk_msg}",
                        confidence_score=0.4
                    )
                    break

            # ═══════════════════════════════════════════════════════════
            # Step 3: Coder - Property Calculation (with Healer retry)
            # Fallback: Direct RDKit if Coder sandbox fails
            # ═══════════════════════════════════════════════════════════
            await self._log_agent_event(session_id, "coder", 3, "started",
                                        "Calculating molecular properties via Coder agent...")

            await self._broadcast_step(session_id, 3, "Coder: Property Calculation...")

            calculated_metrics = {}
            calc_warnings = []
            critical_errors = []

            try:
                # PRIMARY: Coder agent execute() with Healer retry
                coder_result = await self._run_coder_with_healer_retry(state, session_id)

                if coder_result.success:
                    calculated_metrics = coder_result.data.get("metrics", {})
                    state["calculated_metrics"] = calculated_metrics
                    await self._log_agent_event(
                        session_id, "coder", 3, "completed",
                        f"Coder calculated {len(calculated_metrics)} properties",
                        reasoning=coder_result.reasoning,
                        confidence_score=coder_result.confidence_score or 0.9
                    )
                else:
                    # Coder+Healer 실패 → Direct RDKit (STRICT, no estimates)
                    logger.warning(f"[navigator-v2] Coder failed: {coder_result.error}, attempting direct RDKit...")
                    await self._log_agent_event(
                        session_id, "coder", 3, "fallback",
                        f"Coder sandbox failed ({coder_result.error}). Attempting direct RDKit..."
                    )
                    try:
                        calculated_metrics, calc_warnings = await self._calculate_properties_direct(
                            state, session_id
                        )
                        state["calculated_metrics"] = calculated_metrics
                        calc_warnings.append(f"Coder agent failed, used direct RDKit (no estimates)")
                    except (ValueError, ImportError) as direct_err:
                        # FAIL-FAST: RDKit도 실패하면 Critical Error
                        error_msg = str(direct_err)
                        critical_errors.append(error_msg)
                        await self._log_agent_event(
                            session_id, "coder", 3, "critical_error", error_msg
                        )
            except Exception as coder_exc:
                # Coder 예외 → Direct RDKit (STRICT)
                logger.warning(f"[navigator-v2] Coder exception: {coder_exc}, attempting direct RDKit...")
                await self._log_agent_event(
                    session_id, "coder", 3, "fallback",
                    f"Coder exception: {str(coder_exc)[:100]}. Attempting direct RDKit..."
                )
                try:
                    calculated_metrics, calc_warnings = await self._calculate_properties_direct(
                        state, session_id
                    )
                    state["calculated_metrics"] = calculated_metrics
                    calc_warnings.append(f"Coder exception, used direct RDKit (no estimates)")
                except (ValueError, ImportError) as direct_err:
                    error_msg = str(direct_err)
                    critical_errors.append(error_msg)
                    await self._log_agent_event(
                        session_id, "coder", 3, "critical_error", error_msg
                    )

            warnings.extend(calc_warnings)

            # ═══════════════════════════════════════════════════════════
            # Step 4: Auditor - Physical Validation
            # ═══════════════════════════════════════════════════════════
            await self._log_agent_event(session_id, "auditor", 4, "started",
                                        "Validating chemical structure via Auditor agent...")

            await self._broadcast_step(session_id, 4, "Auditor: Physical Validation...")

            physics_verified = False
            validation_flags = {}

            if not calculated_metrics:
                # Step 3에서 물성 연산 자체가 실패한 경우
                critical_errors.append(
                    "물성 검증 불가: Step 3 Property Calculation이 실패하여 "
                    "검증할 데이터가 없습니다."
                )
                await self._log_agent_event(
                    session_id, "auditor", 4, "critical_error",
                    "No calculated metrics to validate — Step 3 failed"
                )
            else:
                try:
                    # PRIMARY: Auditor agent execute()
                    auditor_result = await self.auditor.execute(state)

                    # IMPORTANT: Auditor의 success는 "승인 여부"이지 "에이전트 실행 성공"이 아님.
                    # data에 decision이 있으면 Auditor가 정상 실행된 것.
                    auditor_has_result = bool(auditor_result.data.get("decision"))

                    if auditor_has_result:
                        validation_flags = auditor_result.data.get("chemistry_validation", {})
                        state["validation_flags"] = validation_flags
                        state["constraint_check"] = auditor_result.data.get("constraint_check", {})
                        physics_verified = auditor_result.data.get("decision", {}).get("approved", False)
                        await self._log_agent_event(
                            session_id, "auditor", 4, "completed",
                            f"Auditor validation: {'APPROVED' if physics_verified else 'REJECTED'} "
                            f"(action: {auditor_result.data.get('decision', {}).get('action', 'N/A')})",
                            reasoning=auditor_result.reasoning,
                            confidence_score=auditor_result.confidence_score
                        )
                    elif auditor_result.error:
                        # 실제 에이전트 실행 에러 (Gemini API 실패 등)
                        validation_flags, physics_verified = self._validate_properties_direct(calculated_metrics)
                        warnings.append(
                            f"Auditor agent 실패 ({auditor_result.error}). "
                            f"기본 범위 검증(MW/LogP/HBD/HBA)으로 대체하였습니다."
                        )
                        await self._log_agent_event(
                            session_id, "auditor", 4, "fallback",
                            f"Auditor failed: {auditor_result.error}. Using range-check validation.",
                            confidence_score=0.4
                        )
                    else:
                        # success=False, error=None, data도 없음 → 알 수 없는 상태
                        validation_flags, physics_verified = self._validate_properties_direct(calculated_metrics)
                        warnings.append(
                            "Auditor agent가 결과 없이 종료되었습니다. "
                            "기본 범위 검증(MW/LogP/HBD/HBA)으로 대체하였습니다."
                        )
                        await self._log_agent_event(
                            session_id, "auditor", 4, "fallback",
                            "Auditor returned no decision data. Using range-check validation.",
                            confidence_score=0.4
                        )
                except Exception as audit_err:
                    logger.warning(f"[navigator-v2] Auditor exception: {audit_err}")
                    validation_flags, physics_verified = self._validate_properties_direct(calculated_metrics)
                    warnings.append(
                        f"Auditor agent 예외 발생: {str(audit_err)[:100]}. "
                        f"기본 범위 검증으로 대체하였습니다."
                    )
                    await self._log_agent_event(
                        session_id, "auditor", 4, "fallback",
                        f"Auditor exception: {str(audit_err)[:100]}. Using range-check validation.",
                        confidence_score=0.3
                    )
            
            # ═══════════════════════════════════════════════════════════
            # Step 5: Virtual Trial Simulation
            # ═══════════════════════════════════════════════════════════
            await self._log_agent_event(session_id, "virtual_trial", 5, "started",
                                        "Running PK/PD simulation...")
            
            await self._broadcast_step(session_id, 5, "Clinical Simulation...")
            
            # Virtual trial 계산 (기존 로직 개선)
            virtual_trial = await self._run_virtual_trial_v2(
                state, antibody_candidates, validation_flags
            )
            
            await self._log_agent_event(
                session_id, "virtual_trial", 5, "completed",
                f"Predicted ORR: {virtual_trial.get('predicted_orr') or 0:.1f}%",
                confidence_score=virtual_trial.get('confidence', 0)
            )
            
            # ═══════════════════════════════════════════════════════════
            # 완료 및 결과 조합
            # ═══════════════════════════════════════════════════════════
            execution_time = (datetime.utcnow() - start_time).total_seconds()
            
            # Data Lineage 수집
            data_lineage = self._collect_data_lineage(
                session_id, disease_name, antibody_candidates, 
                self.agent_logs, state
            )
            
            await self._log_agent_event(session_id, "orchestrator", 5, "completed",
                                        f"Navigator completed in {execution_time:.1f}s")
            
            # Virtual trial에 에러가 있으면 critical_errors에 추가
            vt_error = virtual_trial.get("error")
            if vt_error:
                critical_errors.append(vt_error)

            return NavigatorResultV2(
                session_id=session_id,
                disease_name=disease_name,
                antibody_candidates=antibody_candidates,
                golden_combination=self._build_golden_combination(state, candidates, antibody_candidates),
                calculated_metrics=state.get("calculated_metrics", {}),
                physics_verified=physics_verified,
                virtual_trial=virtual_trial,
                agent_logs=self.agent_logs,
                data_lineage=data_lineage,
                execution_time_seconds=execution_time,
                warnings=warnings,
                critical_errors=critical_errors,
            )

        except Exception as e:
            logger.exception(f"[navigator-v2] Pipeline CRITICAL error: {e}")
            await self._log_agent_event(session_id, "orchestrator", state.get("step", 0),
                                        "critical_error", str(e))

            # FAIL-FAST: 파이프라인 전체 실패 시 명확한 에러 + 부분 결과
            return NavigatorResultV2(
                session_id=session_id,
                disease_name=disease_name,
                antibody_candidates=antibody_candidates if 'antibody_candidates' in locals() else [],
                golden_combination={},
                calculated_metrics=state.get("calculated_metrics", {}),
                physics_verified=False,
                virtual_trial={},
                agent_logs=self.agent_logs,
                data_lineage={"error": str(e)},
                execution_time_seconds=(datetime.utcnow() - start_time).total_seconds(),
                warnings=warnings,
                critical_errors=[f"Pipeline Critical Error: {str(e)}"],
            )

    async def _run_coder_with_healer_retry(
        self, 
        state: DesignSessionState, 
        session_id: str,
        max_retries: int = 3
    ) -> AgentOutput:
        """
        FIXED: Coder 실행 with Healer 자동 재시도
        """
        retry_count = 0
        
        while retry_count < max_retries:
            # Coder 실행
            coder_result = await self.coder.execute(state)
            
            if coder_result.success:
                # 성공 시 metrics 저장
                if coder_result.data.get("metrics"):
                    state["calculated_metrics"] = coder_result.data["metrics"]
                return coder_result
            
            # 실패 시 Healer 호출
            retry_count += 1
            if retry_count >= max_retries:
                break
            
            await self._log_agent_event(
                session_id, "healer", 3, "healing",
                f"Fixing Coder error (attempt {retry_count}/{max_retries}): {coder_result.error}"
            )
            
            await self._broadcast_step(session_id, 3, f"Healer: Auto-fixing... ({retry_count})")
            
            # Healer 실행
            state["requires_healing"] = True
            state["healing_attempts"] = retry_count
            state["last_error"] = coder_result.error
            state["last_code"] = coder_result.data.get("code", "")
            
            healer_result = await self.healer.execute(state)
            
            if healer_result.success:
                # Healer가 수정한 코드로 재시도
                state["requires_healing"] = False
                fixed_code = healer_result.data.get("fixed_code")
                if fixed_code:
                    await self._log_agent_event(
                        session_id, "healer", 3, "completed",
                        f"Applied fix: {healer_result.reasoning[:100]}..."
                    )
            else:
                await self._log_agent_event(
                    session_id, "healer", 3, "error",
                    f"Failed to fix: {healer_result.error}"
                )
                break
        
        # 최종 실패 반환
        return AgentOutput(
            success=False,
            data={},
            error=f"Failed after {retry_count} healing attempts"
        )

    # 질환별 대표 항체 매핑 (golden_set 기반)
    DISEASE_ANTIBODY_MAP: Dict[str, List[Dict[str, Any]]] = {
        "breast cancer": [
            {"name": "Trastuzumab", "target_protein": "HER2", "clinical_score": 0.95},
            {"name": "Sacituzumab", "target_protein": "TROP2", "clinical_score": 0.88},
            {"name": "Patritumab", "target_protein": "HER3", "clinical_score": 0.72},
        ],
        "lung cancer": [
            {"name": "Trastuzumab", "target_protein": "HER2", "clinical_score": 0.82},
            {"name": "Datopotamab", "target_protein": "TROP2", "clinical_score": 0.80},
            {"name": "Telisotuzumab", "target_protein": "c-Met", "clinical_score": 0.70},
        ],
        "urothelial cancer": [
            {"name": "Enfortumab", "target_protein": "Nectin-4", "clinical_score": 0.90},
            {"name": "Sacituzumab", "target_protein": "TROP2", "clinical_score": 0.85},
        ],
        "lymphoma": [
            {"name": "Brentuximab", "target_protein": "CD30", "clinical_score": 0.92},
            {"name": "Polatuzumab", "target_protein": "CD79b", "clinical_score": 0.85},
            {"name": "Loncastuximab", "target_protein": "CD19", "clinical_score": 0.80},
        ],
        "aml": [
            {"name": "Gemtuzumab", "target_protein": "CD33", "clinical_score": 0.78},
        ],
        "gastric cancer": [
            {"name": "Trastuzumab", "target_protein": "HER2", "clinical_score": 0.88},
            {"name": "Zolbetuximab", "target_protein": "CLDN18.2", "clinical_score": 0.75},
        ],
        "ovarian cancer": [
            {"name": "Mirvetuximab", "target_protein": "FRa", "clinical_score": 0.86},
        ],
        "melanoma": [
            {"name": "Glembatumumab", "target_protein": "gpNMB", "clinical_score": 0.55},
        ],
    }

    # gpNMB 등 임상 실패 타겟 목록
    HIGH_RISK_TARGETS = {
        "gpNMB": "METRIC trial failed (Phase II): Glembatumumab vedotin showed no PFS benefit vs. capecitabine in gpNMB+ TNBC (HR=0.95, p=0.83)",
        "EGFR": "Depatuxizumab mafodotin (ABT-414) failed Phase III INTELLANCE-1 in GBM",
        "MUC16": "Sofituzumab vedotin discontinued in ovarian cancer (limited efficacy)",
    }

    # PBD (Pyrrolobenzodiazepine) 등 고독성 payload 클래스
    HIGH_TOXICITY_PAYLOADS = {
        "PBD": "PBD dimers carry severe toxicity risk: Vadastuximab talirine (SGN-CD33A) discontinued after fatal hepatotoxicity in CASCADE trial. Loncastuximab tesirine has narrow therapeutic window.",
        "PBD dimer": "PBD dimers carry severe toxicity risk: Vadastuximab talirine (SGN-CD33A) discontinued after fatal hepatotoxicity in CASCADE trial.",
        "Calicheamicin": "Calicheamicin-based ADCs (Gemtuzumab ozogamicin) have hepatic veno-occlusive disease risk. Require careful dose management.",
    }

    def _extract_antibodies_from_librarian(
        self,
        result: AgentOutput,
        disease_name: str
    ) -> List[Dict[str, Any]]:
        """Librarian 결과에서 항체 목록 추출 (질환별 매핑)"""
        candidates = []

        # Golden Set 참조에서 항체 정보 추출
        gs_refs = result.data.get("golden_set_references", [])
        for ref in gs_refs[:3]:
            candidates.append({
                "id": ref.get("id", str(uuid.uuid4())[:8]),
                "name": ref.get("name", "Unknown"),
                "target_protein": ref.get("target", disease_name),
                "clinical_score": ref.get("clinical_score", 0.85),
                "match_confidence": result.confidence_score or 0.8,
                "source": "librarian_golden_set"
            })

        # 부족하면 질환별 매핑에서 보충
        if len(candidates) < 3:
            disease_lower = disease_name.lower().strip()
            fallback_abs = []

            # 부분 매칭 (e.g. "triple negative breast cancer" → "breast cancer")
            for disease_key, abs_list in self.DISEASE_ANTIBODY_MAP.items():
                if disease_key in disease_lower or disease_lower in disease_key:
                    fallback_abs = abs_list
                    break

            existing_names = {c["name"] for c in candidates}
            for ab in fallback_abs:
                if ab["name"] not in existing_names and len(candidates) < 3:
                    candidates.append({
                        "id": f"{ab['name'].lower()}-fallback",
                        "name": ab["name"],
                        "target_protein": ab["target_protein"],
                        "clinical_score": ab["clinical_score"],
                        "match_confidence": 0.70,
                        "source": "disease_mapping_fallback"
                    })

        return candidates[:3]

    def _build_golden_combination(
        self,
        state: DesignSessionState,
        candidates: List[CandidateStructure],
        antibody_candidates: List[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        최종 조합 구성 — golden_set_library 실제 데이터 필수

        FAIL-FAST: golden_set에서 실제 레퍼런스를 찾지 못하면
        데이터 출처를 명시하고, 하드코딩 기본값 사용 금지.
        """
        target = state.get("target_antigen", "")
        antibody_name = None
        linker_type = None
        payload_class = None
        dar = None
        historical_orr = None
        data_source = "unknown"

        # 1순위: golden_set_library에서 FDA 승인 레퍼런스 검색
        try:
            gs_ref = self.supabase.table("golden_set_library").select(
                "name, target_1, linker_type, dar, orr_pct, description, payload_class"
            ).eq("status", "approved").eq("target_1", target).limit(1).execute()

            if gs_ref.data:
                ref = gs_ref.data[0]
                antibody_name = ref.get("name")
                linker_type = ref.get("linker_type")
                dar = ref.get("dar")
                historical_orr = ref.get("orr_pct")
                payload_class = ref.get("payload_class")
                data_source = f"golden_set_library (approved, target={target})"
                logger.info(f"[navigator-v2] Golden combination from DB: {antibody_name}, {linker_type}, DAR={dar}")
        except Exception as e:
            logger.error(f"[navigator-v2] Golden set lookup FAILED: {e}")

        # 2순위: antibody_candidates에서 항체명 보충
        if not antibody_name and antibody_candidates:
            top = antibody_candidates[0]
            antibody_name = top.get("name")
            data_source = f"librarian_candidates (target={target})"

        # Alchemist 결과에서 SMILES/payload 보충
        linker_smiles = ""
        payload_smiles = ""
        if candidates:
            top_candidate = candidates[0]
            linker_smiles = top_candidate.get("smiles", "")
            payload_smiles = top_candidate.get("payload_smiles", linker_smiles)
            if not payload_class and top_candidate.get("payload_class"):
                payload_class = top_candidate["payload_class"]

        # FAIL-FAST: 필수 필드 누락 시 명시적 경고 (가짜값 대입 안함)
        missing_fields = []
        if not antibody_name:
            missing_fields.append("antibody")
        if not linker_type:
            missing_fields.append("linker_type")
        if not payload_class:
            missing_fields.append("payload_class")
        if not dar:
            missing_fields.append("DAR")

        if missing_fields:
            logger.warning(
                f"[navigator-v2] INCOMPLETE golden combination for target '{target}': "
                f"missing {missing_fields}. No fake defaults applied."
            )

        return {
            "antibody": antibody_name,
            "linker": {
                "type": linker_type,
                "smiles": linker_smiles
            },
            "payload": {
                "class": payload_class,
                "smiles": payload_smiles
            },
            "dar": dar,
            "historical_orr": historical_orr,
            "data_source": data_source,
            "missing_fields": missing_fields if missing_fields else None,
        }

    async def _calculate_properties_direct(
        self,
        state: DesignSessionState,
        session_id: str
    ) -> tuple:
        """
        RDKit 직접 계산 (Coder+Healer 실패 시 최후의 보루)

        FAIL-FAST 원칙:
        - SMILES 없음 → CRITICAL ERROR (추정값 절대 사용 금지)
        - RDKit 미설치 → CRITICAL ERROR (환경 문제 즉시 노출)
        - Invalid SMILES → CRITICAL ERROR (잘못된 구조 즉시 노출)
        """
        warnings = []
        metrics = {}

        smiles = state.get("current_smiles", "")
        if not smiles:
            candidates = state.get("candidates", [])
            if candidates:
                smiles = candidates[0].get("smiles", "")

        if not smiles:
            raise ValueError(
                "CRITICAL: No SMILES structure available. "
                "Alchemist가 유효한 분자 구조를 생성하지 못했습니다. "
                "물성 연산을 수행할 수 없습니다."
            )

        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors
        except ImportError:
            raise ImportError(
                "CRITICAL: RDKit 물리 엔진이 설치되지 않았습니다. "
                "Cloud Run 환경에 rdkit>=2024.3.6이 누락되어 정확한 물성 연산이 불가능합니다. "
                "추정값 사용은 허용되지 않습니다. 즉시 환경을 점검하십시오."
            )

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(
                f"CRITICAL: Invalid SMILES 구조입니다 — '{smiles[:80]}'. "
                f"RDKit이 파싱할 수 없는 잘못된 분자식입니다. "
                f"Alchemist가 생성한 SMILES를 검증하십시오."
            )

        metrics = {
            "mw": round(Descriptors.MolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "tpsa": round(Descriptors.TPSA(mol), 2),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "num_rings": rdMolDescriptors.CalcNumRings(mol),
            "smiles": smiles,
            "computed_by": "rdkit_direct",
        }

        # SA Score (선택적 — 이것만 실패해도 전체는 성공)
        try:
            from rdkit.Chem import RDConfig
            import os, sys
            sa_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
            if os.path.isdir(sa_path):
                sys.path.insert(0, sa_path)
                import sascorer
                metrics["sa_score"] = round(sascorer.calculateScore(mol), 2)
        except Exception:
            warnings.append("SA Score calculation unavailable (non-critical)")

        return metrics, warnings

    def _validate_properties_direct(self, metrics: Dict) -> tuple:
        """직접 물성 검증 (Auditor 실패 시 fallback)"""
        flags = {}
        mw = metrics.get("mw", 0)
        logp = metrics.get("logp", 0)
        hbd = metrics.get("hbd", 0)
        hba = metrics.get("hba", 0)

        # ADC payload는 Lipinski rule 범위 밖일 수 있지만 기본 체크
        flags["mw_range"] = 200 < mw < 2000 if mw else True
        flags["logp_range"] = -2 < logp < 8 if logp else True
        flags["hbd_check"] = hbd <= 10 if hbd else True
        flags["hba_check"] = hba <= 20 if hba else True
        flags["lipinski_pass"] = all([
            mw < 1000, logp < 6, hbd <= 5, hba <= 10
        ]) if mw else False

        # ADC payload는 대부분 Lipinski 위반이므로 다른 기준 적용
        verified = flags.get("mw_range", True) and flags.get("logp_range", True)

        return flags, verified

    async def _run_virtual_trial_v2(
        self,
        state: DesignSessionState,
        antibody_candidates: List[Dict],
        validation_flags: Dict
    ) -> Dict[str, Any]:
        """
        Virtual Trial 시뮬레이션 V2 — golden_set 실제 임상 데이터 기반

        FAIL-FAST 원칙:
        - golden_set에서 실제 ORR/PFS/OS 데이터를 조회
        - 데이터가 없으면 '데이터 부족으로 예측 불가'라고 정직하게 보고
        - 하드코딩 ORR=45, PFS=6 등 절대 사용 금지
        """
        target = state.get("target_antigen", "")
        disease = state.get("target_indication", "")

        # golden_set_library에서 동일 타겟의 실제 임상 데이터 조회
        ref_orr = None
        ref_pfs = None
        ref_os = None
        ref_name = None
        data_source = None

        try:
            gs_refs = self.supabase.table("golden_set_library").select(
                "name, target_1, orr_pct, pfs_months, os_months, clinical_status, indication"
            ).eq("target_1", target).not_.is_("orr_pct", "null").order(
                "orr_pct", desc=True
            ).limit(3).execute()

            if gs_refs.data:
                # 가장 높은 ORR 레퍼런스 사용
                best_ref = gs_refs.data[0]
                ref_orr = best_ref.get("orr_pct")
                ref_pfs = best_ref.get("pfs_months")
                ref_os = best_ref.get("os_months")
                ref_name = best_ref.get("name")
                data_source = f"golden_set_library: {ref_name} (target={target})"
                logger.info(f"[navigator-v2] Virtual trial ref: {ref_name} ORR={ref_orr}%, PFS={ref_pfs}m, OS={ref_os}m")
        except Exception as e:
            logger.error(f"[navigator-v2] Golden set clinical data query FAILED: {e}")

        # FAIL-FAST: 실제 레퍼런스 데이터가 없으면 정직하게 보고
        if ref_orr is None:
            logger.warning(f"[navigator-v2] No clinical reference data for target '{target}' — prediction unavailable")
            return {
                "predicted_orr": None,
                "predicted_pfs_months": None,
                "predicted_os_months": None,
                "confidence": 0.0,
                "pk_data": [],
                "tumor_data": [],
                "data_source": None,
                "reference_drug": None,
                "error": (
                    f"임상 예측 데이터 부족: target '{target}'에 대한 "
                    f"Golden Set 레퍼런스 ORR/PFS/OS 데이터가 존재하지 않습니다. "
                    f"가상 임상 시뮬레이션을 수행할 수 없습니다."
                ),
            }

        # Validation 결과에 따른 미세 조정 (레퍼런스 기반)
        adjustment = 1.0
        if validation_flags.get("lipinski_pass"):
            adjustment *= 1.05  # 물성 통과 시 소폭 보너스
        if validation_flags.get("pains_detected"):
            adjustment *= 0.85  # PAINS 검출 시 감점

        # 신뢰도 계산
        confidence = 0.5  # 기본 (레퍼런스 존재)
        if ref_pfs is not None:
            confidence += 0.15
        if ref_os is not None:
            confidence += 0.15
        if len(gs_refs.data) >= 2:
            confidence += 0.1  # 다수 레퍼런스 존재

        predicted_orr = round(min(ref_orr * adjustment, 95), 1)
        predicted_pfs = round(ref_pfs * adjustment, 1) if ref_pfs else None
        predicted_os = round(ref_os * adjustment, 1) if ref_os else None

        # PK/Tumor 시뮬레이션은 레퍼런스 기반으로 생성
        half_life_hours = 96 if ref_pfs and ref_pfs > 6 else 72  # 레퍼런스 기반 추정
        pk_data = self._generate_pk_data(half_life_hours=half_life_hours)
        tumor_data = self._generate_tumor_data(predicted_orr)

        return {
            "predicted_orr": predicted_orr,
            "predicted_pfs_months": predicted_pfs,
            "predicted_os_months": predicted_os,
            "confidence": round(confidence, 2),
            "pk_data": pk_data,
            "tumor_data": tumor_data,
            "data_source": data_source,
            "reference_drug": ref_name,
            "reference_count": len(gs_refs.data) if gs_refs.data else 0,
        }

    def _generate_pk_data(self, half_life_hours: float = 96) -> List[Dict]:
        """PK 데이터 생성 (레퍼런스 반감기 기반)"""
        # half_life를 일 단위 decay rate로 변환
        decay_per_day = 0.5 ** (24 / half_life_hours)
        return [
            {
                "time_hours": i * 24,
                "concentration": round(100 * (decay_per_day ** i), 2),
                "free_payload": round(5 * (decay_per_day ** i), 3),
            }
            for i in range(15)
        ]

    def _generate_tumor_data(self, orr: float) -> List[Dict]:
        """종양 데이터 생성"""
        tgi = min(orr * 1.1, 95)
        return [
            {
                "day": i * 3,
                "control": 100 * (2 ** (i * 3 / 5)),
                "treated": 100 * (2 ** (i * 3 / 5)) * (1 - tgi / 100 * (1 - 0.5 ** (i / 7)))
            }
            for i in range(15)
        ]

    def _collect_data_lineage(
        self,
        session_id: str,
        disease_name: str,
        antibody_candidates: List[Dict],
        agent_logs: List[NavigatorAgentLog],
        state: DesignSessionState
    ) -> Dict[str, Any]:
        """데이터 출처 추적 정보 수집"""
        return {
            "pipeline_version": "2.0-agent-integrated",
            "session_id": session_id,
            "disease_input": disease_name,
            "execution_timestamp": datetime.utcnow().isoformat(),
            "agents_invoked": ["librarian", "alchemist", "coder", "healer", "auditor"],
            "agent_execution_order": [log.agent_name for log in agent_logs],
            "data_sources": {
                "golden_set_references": state.get("golden_set_references", []),
                "pmid_references": list(set([
                    pmid for log in agent_logs 
                    for pmid in (log.pmid_references or [])
                ])),
                "reagent_references": state.get("scaffold_mappings", []),
            },
            "agent_logs": [
                {
                    "agent": log.agent_name,
                    "step": log.step,
                    "status": log.status,
                    "data_source": log.data_source,
                    "pmid_refs": log.pmid_references,
                    "confidence": log.confidence_score
                }
                for log in agent_logs
            ]
        }

    async def _log_agent_event(
        self,
        session_id: str,
        agent_name: str,
        step: int,
        status: str,
        message: str,
        reasoning: Optional[str] = None,
        data_source: Optional[str] = None,
        pmid_references: Optional[List[str]] = None,
        confidence_score: Optional[float] = None
    ):
        """에이전트 이벤트 로깅 및 스트리밍"""
        log = NavigatorAgentLog(
            timestamp=datetime.utcnow().isoformat(),
            agent_name=agent_name,
            step=step,
            status=status,
            message=message,
            reasoning=reasoning,
            data_source=data_source,
            pmid_references=pmid_references or [],
            confidence_score=confidence_score
        )
        
        self.agent_logs.append(log)
        
        # WebSocket 스트리밍
        try:
            await websocket_hub.stream_agent_log(
                session_id=session_id,
                level="info" if status in ["started", "completed"] else "warning",
                message=message,
                emoji=self._get_agent_emoji(agent_name),
                agent_name=agent_name,
                metadata={
                    "step": step,
                    "status": status,
                    "reasoning": reasoning,
                    "data_source": data_source
                }
            )
        except Exception as e:
            logger.warning(f"[navigator-v2] WebSocket broadcast error: {e}")
        
        # DB에도 저장 (모든 에이전트 로그 기록)
        try:
            self.supabase.table("agent_execution_logs").insert({
                "session_id": session_id,
                "agent_name": agent_name,
                "status": status,
                "reasoning": reasoning,
                "decision_summary": message,
                "confidence_score": confidence_score,
                "referenced_pmids": pmid_references,
                "created_at": datetime.utcnow().isoformat()
            }).execute()
        except Exception as e:
            logger.warning(f"[navigator-v2] DB log error for {agent_name}: {e}")

    def _get_agent_emoji(self, agent_name: str) -> str:
        """에이전트별 이모지"""
        emojis = {
            "librarian": "📚",
            "alchemist": "⚗️",
            "coder": "💻",
            "healer": "🔧",
            "auditor": "🔍",
            "virtual_trial": "🏥",
            "orchestrator": "🎭"
        }
        return emojis.get(agent_name, "🤖")

    async def _broadcast_step(self, session_id: str, step: int, message: str):
        """진행 상황 브로드캐스트"""
        try:
            await websocket_hub.stream_progress(
                session_id=session_id,
                percentage=step * 20,
                operation=message
            )
        except Exception as e:
            logger.warning(f"[navigator-v2] Progress broadcast error: {e}")


# ============================================================================
# Public API
# ============================================================================

_navigator_v2: Optional[NavigatorOrchestratorV2] = None


def get_navigator_orchestrator_v2() -> NavigatorOrchestratorV2:
    """Navigator V2 싱글톤"""
    global _navigator_v2
    if _navigator_v2 is None:
        _navigator_v2 = NavigatorOrchestratorV2()
    return _navigator_v2


async def run_one_click_navigator_v2(
    disease_name: str,
    session_id: Optional[str] = None,
    user_id: Optional[str] = None,
    target_protein: Optional[str] = None,
    selected_antibody_id: Optional[str] = None
) -> NavigatorResultV2:
    """
    One-Click ADC Navigator V2 실행

    실제 에이전트 호출 + DesignSessionState 공유 + 타겟 선택
    """
    orchestrator = get_navigator_orchestrator_v2()
    return await orchestrator.run_navigator_v2(
        disease_name=disease_name,
        session_id=session_id,
        user_id=user_id,
        target_protein=target_protein,
        selected_antibody_id=selected_antibody_id
    )
