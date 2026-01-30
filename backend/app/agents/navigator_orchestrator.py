"""
One-Click ADC Navigator Orchestrator v3.0
AstraForge v2.2 - Real Agent Integration

기존 Facade 패턴을 제거하고, DesignOrchestrator에서 검증된
6인 AI 에이전트 프레임워크(execute() 메서드)를 직접 호출합니다.

파이프라인:
- Step 1: Target Discovery (Librarian.execute → target/antibody 발견)
- Step 2: Assembly (Alchemist.execute → Golden Combination SMILES 생성)
- Step 3: Calculation (Coder.execute + Healer → Sandbox 물성 계산)
- Step 4: Validation (Auditor.execute → PAINS/Lipinski/독성 검증)
- Step 5: Clinical Prediction (Gemini → 실제 golden_set 기반 임상 예측)
"""

import logging
import json
import uuid
import time
import math
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum

from app.core.supabase import get_supabase_client
from app.core.websocket_hub import websocket_hub
from app.core.gemini import get_gemini_model
from app.agents.design_state import (
    DesignSessionState,
    create_initial_state,
)
from app.agents.base_agent import AgentOutput

logger = logging.getLogger(__name__)


# ============================================================================
# Data Classes (프론트엔드 API 계약 유지)
# ============================================================================

class NavigatorStep(Enum):
    INIT = 0
    TARGET_MATCH = 1
    GOLDEN_COMBINATION = 2
    PROPERTY_CALCULATION = 3
    PHYSICAL_VALIDATION = 4
    VIRTUAL_TRIAL = 5
    COMPLETE = 6


@dataclass
class AntibodyCandidate:
    antibody_id: str
    name: str
    target_protein: str
    isotype: Optional[str] = None
    related_diseases: Optional[str] = None
    full_spec: Optional[Dict[str, Any]] = None
    clinical_score: float = 0.0
    match_confidence: float = 0.0


@dataclass
class LinkerSpec:
    id: str
    smiles: str
    type: str
    cleavable: bool = True
    properties: Optional[Dict[str, Any]] = None


@dataclass
class PayloadSpec:
    id: str
    smiles: str
    class_name: str
    mechanism: Optional[str] = None
    properties: Optional[Dict[str, Any]] = None


@dataclass
class GoldenCombination:
    antibody: AntibodyCandidate
    linker: LinkerSpec
    payload: PayloadSpec
    dar: int = 4
    historical_performance: Optional[Dict[str, Any]] = None
    confidence_score: float = 0.0


@dataclass
class VirtualTrialResult:
    predicted_orr: float
    predicted_pfs_months: float
    predicted_os_months: float
    pk_data: List[Dict[str, float]] = field(default_factory=list)
    tumor_data: List[Dict[str, float]] = field(default_factory=list)
    patient_population: str = ""
    confidence: float = 0.7


@dataclass
class NavigatorResult:
    session_id: str
    disease_name: str
    antibody_candidates: List[AntibodyCandidate]
    golden_combination: GoldenCombination
    calculated_metrics: Dict[str, Any]
    physical_validations: List[Dict[str, Any]]
    physics_verified: bool
    virtual_trial: VirtualTrialResult
    digital_lineage: Dict[str, Any]
    combined_smiles: str = ""
    execution_time_seconds: float = 0.0


@dataclass
class NavigatorState:
    session_id: str
    user_id: Optional[str]
    disease_name: str
    step: NavigatorStep = NavigatorStep.INIT
    total_steps: int = 5
    antibody_candidates: List[AntibodyCandidate] = field(default_factory=list)
    target_protein: str = ""
    golden_combination: Optional[GoldenCombination] = None
    combined_smiles: str = ""
    calculated_metrics: Dict[str, Any] = field(default_factory=dict)
    simulation_results: Dict[str, Any] = field(default_factory=dict)
    physical_validations: List[Dict[str, Any]] = field(default_factory=list)
    physics_verified: bool = False
    virtual_trial: Optional[VirtualTrialResult] = None
    started_at: datetime = field(default_factory=datetime.utcnow)
    errors: List[str] = field(default_factory=list)


# ============================================================================
# Disease Name Mapping
# ============================================================================

DISEASE_NAME_MAP = {
    "유방암": "Breast Cancer",
    "폐암": "Lung Cancer",
    "위암": "Gastric Cancer",
    "대장암": "Colorectal Cancer",
    "간암": "Liver Cancer",
    "췌장암": "Pancreatic Cancer",
    "방광암": "Bladder Cancer",
    "난소암": "Ovarian Cancer",
    "자궁경부암": "Cervical Cancer",
    "전립선암": "Prostate Cancer",
    "백혈병": "Leukemia",
    "림프종": "Lymphoma",
    "다발성골수종": "Multiple Myeloma",
    "흑색종": "Melanoma",
    "두경부암": "Head and Neck Cancer",
    "신장암": "Renal Cell Carcinoma",
    "갑상선암": "Thyroid Cancer",
}


# ============================================================================
# Navigator Orchestrator v3 - Real Agent Integration
# ============================================================================

class NavigatorOrchestrator:
    """
    One-Click ADC Navigator v3.0 - Real Agent Integration

    Facade 패턴을 완전히 제거하고 실제 에이전트의 execute()를 호출합니다.
    - Step 1: LibrarianAgent.execute() + find_antibodies_by_disease()
    - Step 2: AlchemistAgent.execute() → Gemini 기반 SMILES 생성
    - Step 3: CoderAgent.execute() → Sandbox 물성 계산
    - Step 4: AuditorAgent.execute() → PAINS/Lipinski/독성 검증
    - Step 5: Gemini 기반 Virtual Trial (golden_set 레퍼런스)
    """

    def __init__(self):
        self.supabase = get_supabase_client()

        # 실제 에이전트 인스턴스 (lazy init)
        self._librarian = None
        self._alchemist = None
        self._coder = None
        self._healer = None
        self._auditor = None

    @property
    def librarian(self):
        if self._librarian is None:
            from app.agents.librarian import LibrarianAgent
            self._librarian = LibrarianAgent()
        return self._librarian

    @property
    def alchemist(self):
        if self._alchemist is None:
            from app.agents.alchemist import AlchemistAgent
            self._alchemist = AlchemistAgent()
        return self._alchemist

    @property
    def coder(self):
        if self._coder is None:
            from app.agents.coder import CoderAgent
            self._coder = CoderAgent()
        return self._coder

    @property
    def healer(self):
        if self._healer is None:
            from app.agents.healer import HealerAgent
            self._healer = HealerAgent()
        return self._healer

    @property
    def auditor(self):
        if self._auditor is None:
            from app.agents.auditor import AuditorAgent
            self._auditor = AuditorAgent()
        return self._auditor

    # =====================================================================
    # Main Pipeline
    # =====================================================================

    async def run_one_click_pipeline(
        self,
        disease_name: str,
        session_id: Optional[str] = None,
        user_id: Optional[str] = None
    ) -> NavigatorResult:
        """
        질환명 기반 원클릭 ADC 설계 파이프라인 (v3 - Real Agents)
        """
        pipeline_start = time.time()

        if not session_id:
            session_id = str(uuid.uuid4())

        state = NavigatorState(
            session_id=session_id,
            user_id=user_id,
            disease_name=disease_name
        )

        # Normalize disease name
        normalized_disease = DISEASE_NAME_MAP.get(disease_name, disease_name)
        if normalized_disease != disease_name:
            logger.info(f"[navigator] Disease normalized: '{disease_name}' → '{normalized_disease}'")

        try:
            await self._create_db_session(state)

            # ══════════════════════════════════════════════════
            # Step 1: Target & Antibody Discovery (Librarian)
            # ══════════════════════════════════════════════════
            state.step = NavigatorStep.TARGET_MATCH
            await self._broadcast_step(state, 1, "Librarian 에이전트가 타겟 및 항체를 검색 중...")
            step1_start = time.time()

            # 1a. Librarian으로 질환에 맞는 항체 검색 (실제 execute)
            antibody_candidates = await self._step1_target_discovery(
                session_id, normalized_disease
            )

            if not antibody_candidates:
                raise ValueError(f"질환 '{normalized_disease}'에 대한 항체 후보를 찾을 수 없습니다.")

            state.antibody_candidates = antibody_candidates
            state.target_protein = antibody_candidates[0].target_protein

            step1_time = time.time() - step1_start
            unique_targets = list(set(ab.target_protein for ab in antibody_candidates))
            await self._update_db_step(state, 1, {
                "antibody_candidates": [
                    {
                        "antibody_id": ab.antibody_id,
                        "name": ab.name,
                        "target_protein": ab.target_protein,
                        "clinical_score": ab.clinical_score,
                        "match_confidence": ab.match_confidence
                    } for ab in antibody_candidates
                ],
                "primary_target": state.target_protein
            })
            await self._broadcast_step(
                state, 1,
                f"✅ Step 1 완료 ({step1_time:.1f}s): "
                f"타겟 {', '.join(unique_targets[:3])}, 항체 {len(antibody_candidates)}개 발견"
            )

            # ══════════════════════════════════════════════════
            # Step 2: Golden Combination Assembly (Alchemist)
            # ══════════════════════════════════════════════════
            state.step = NavigatorStep.GOLDEN_COMBINATION
            await self._broadcast_step(state, 2, "Alchemist 에이전트가 최적 조합을 설계 중...")
            step2_start = time.time()

            golden_combo, combined_smiles = await self._step2_assembly(
                session_id, state.target_protein, normalized_disease, antibody_candidates
            )

            state.golden_combination = golden_combo
            state.combined_smiles = combined_smiles

            step2_time = time.time() - step2_start
            await self._update_db_step(state, 2, {
                "golden_combination": {
                    "antibody": golden_combo.antibody.name,
                    "linker_type": golden_combo.linker.type,
                    "payload_class": golden_combo.payload.class_name,
                    "dar": golden_combo.dar,
                },
                "combined_smiles": combined_smiles
            })
            await self._broadcast_step(
                state, 2,
                f"✅ Step 2 완료 ({step2_time:.1f}s): "
                f"{golden_combo.antibody.name} + {golden_combo.payload.class_name} (DAR={golden_combo.dar})"
            )

            # ══════════════════════════════════════════════════
            # Step 3: Property Calculation (Coder + Healer)
            # ══════════════════════════════════════════════════
            state.step = NavigatorStep.PROPERTY_CALCULATION
            await self._broadcast_step(state, 3, "Coder 에이전트가 Sandbox에서 물성을 계산 중...")
            step3_start = time.time()

            calculated_metrics = await self._step3_property_calculation(
                session_id, combined_smiles, state.target_protein, normalized_disease
            )

            state.calculated_metrics = calculated_metrics

            step3_time = time.time() - step3_start
            mw = calculated_metrics.get("mw", "N/A")
            logp = calculated_metrics.get("logp", "N/A")
            await self._update_db_step(state, 3, {"calculated_metrics": calculated_metrics})
            await self._broadcast_step(
                state, 3,
                f"✅ Step 3 완료 ({step3_time:.1f}s): MW={mw}, LogP={logp}"
            )

            # ══════════════════════════════════════════════════
            # Step 4: Physical Validation (Auditor)
            # ══════════════════════════════════════════════════
            state.step = NavigatorStep.PHYSICAL_VALIDATION
            await self._broadcast_step(state, 4, "Auditor 에이전트가 구조 검증 및 실패 사례 체크 중...")
            step4_start = time.time()

            physics_verified, validations, failure_warnings = await self._step4_validation(
                session_id, combined_smiles, state.target_protein,
                normalized_disease, calculated_metrics
            )

            state.physics_verified = physics_verified
            state.physical_validations = validations

            step4_time = time.time() - step4_start
            verification_label = "Physics Verified ✅" if physics_verified else "검증 경고 ⚠️"
            await self._update_db_step(state, 4, {
                "physics_verified": physics_verified,
                "physical_validations": validations,
                "failure_warnings": failure_warnings
            })
            await self._broadcast_step(
                state, 4,
                f"✅ Step 4 완료 ({step4_time:.1f}s): {verification_label}"
                + (f", 경고: {', '.join(failure_warnings[:2])}" if failure_warnings else "")
            )

            # ══════════════════════════════════════════════════
            # Step 5: Virtual Clinical Trial (Gemini + Golden Set)
            # ══════════════════════════════════════════════════
            state.step = NavigatorStep.VIRTUAL_TRIAL
            await self._broadcast_step(state, 5, "실제 임상 데이터 기반 가상 임상 시뮬레이션 중...")
            step5_start = time.time()

            virtual_trial = await self._step5_virtual_trial(
                session_id, golden_combo, calculated_metrics,
                normalized_disease, state.target_protein
            )

            state.virtual_trial = virtual_trial

            step5_time = time.time() - step5_start
            await self._update_db_step(state, 5, {
                "virtual_trial": {
                    "predicted_orr": virtual_trial.predicted_orr,
                    "predicted_pfs_months": virtual_trial.predicted_pfs_months,
                    "predicted_os_months": virtual_trial.predicted_os_months,
                    "pk_data": virtual_trial.pk_data,
                    "tumor_data": virtual_trial.tumor_data,
                    "confidence": virtual_trial.confidence
                }
            })
            await self._broadcast_step(
                state, 5,
                f"✅ Step 5 완료 ({step5_time:.1f}s): "
                f"예측 ORR={virtual_trial.predicted_orr:.1f}%, "
                f"PFS={virtual_trial.predicted_pfs_months:.1f}mo"
            )

            # ══════════════════════════════════════════════════
            # 완료
            # ══════════════════════════════════════════════════
            total_time = time.time() - pipeline_start

            # Digital lineage
            lineage = {
                "version": "v3.0.0-real-agents",
                "pipeline": "navigator_orchestrator",
                "agents_used": ["librarian", "alchemist", "coder", "healer", "auditor"],
                "agent_execution_mode": "execute()",
                "step_times": {
                    "step1_target_discovery": round(step1_time, 2),
                    "step2_assembly": round(step2_time, 2),
                    "step3_calculation": round(step3_time, 2),
                    "step4_validation": round(step4_time, 2),
                    "step5_virtual_trial": round(step5_time, 2),
                },
                "total_execution_seconds": round(total_time, 2),
                "disease_input": disease_name,
                "disease_normalized": normalized_disease,
                "timestamp": datetime.utcnow().isoformat()
            }

            result = NavigatorResult(
                session_id=session_id,
                disease_name=disease_name,
                antibody_candidates=state.antibody_candidates,
                golden_combination=state.golden_combination,
                calculated_metrics=state.calculated_metrics,
                physical_validations=state.physical_validations,
                physics_verified=state.physics_verified,
                virtual_trial=state.virtual_trial,
                digital_lineage=lineage,
                combined_smiles=state.combined_smiles,
                execution_time_seconds=round(total_time, 2)
            )

            logger.info(
                f"[navigator] Pipeline complete in {total_time:.1f}s "
                f"for '{disease_name}' → target={state.target_protein}, "
                f"ORR={virtual_trial.predicted_orr:.1f}%"
            )

            return result

        except Exception as e:
            logger.exception(f"[navigator] Pipeline error: {e}")
            state.errors.append(str(e))
            await self._update_db_error(state, str(e))
            raise

    # =====================================================================
    # Step 1: Target Discovery via Librarian Agent
    # =====================================================================

    async def _step1_target_discovery(
        self, session_id: str, disease_name: str
    ) -> List[AntibodyCandidate]:
        """
        Step 1: Librarian 에이전트의 find_antibodies_by_disease()를 사용하여
        실제 벡터 검색 + golden_set 매칭으로 항체 후보를 찾습니다.
        """
        logger.info(f"[navigator:step1] Calling Librarian for disease: {disease_name}")

        # Librarian의 전용 Navigator 메서드 호출 (실제 RAG + 벡터 검색)
        raw_antibodies = await self.librarian.find_antibodies_by_disease(
            disease_name=disease_name,
            top_k=5
        )

        if not raw_antibodies:
            logger.warning(f"[navigator:step1] Librarian returned no results, trying golden_set fallback")
            raw_antibodies = await self._golden_set_antibody_fallback(disease_name)

        if not raw_antibodies:
            logger.warning(f"[navigator:step1] Golden set fallback also empty, using Gemini")
            raw_antibodies = await self._gemini_target_suggestion(disease_name)

        # Convert to AntibodyCandidate dataclass
        candidates = []
        for ab in raw_antibodies:
            candidate = AntibodyCandidate(
                antibody_id=str(ab.get("id", str(uuid.uuid4()))),
                name=ab.get("name") or ab.get("product_name", "Unknown"),
                target_protein=ab.get("target_protein") or ab.get("target_normalized") or ab.get("target_1", "Unknown"),
                isotype=ab.get("isotype"),
                related_diseases=ab.get("related_disease") or ab.get("category"),
                clinical_score=float(ab.get("clinical_score", 0) or 0),
                match_confidence=float(ab.get("combined_score") or ab.get("similarity", 0) or 0)
            )
            candidates.append(candidate)

        # Sort by clinical score
        candidates.sort(key=lambda x: x.clinical_score, reverse=True)

        # Check for failure cases and add warnings
        for candidate in candidates:
            await self._check_failure_cases(session_id, candidate.target_protein)

        logger.info(f"[navigator:step1] Found {len(candidates)} candidates: "
                    f"{[c.target_protein for c in candidates[:3]]}")

        return candidates[:5]

    async def _golden_set_antibody_fallback(self, disease_name: str) -> List[Dict]:
        """Golden set에서 직접 항체 검색 (Librarian 벡터검색 실패 시)"""
        try:
            search_terms = self._get_search_terms(disease_name)
            all_results = []

            for term in search_terms:
                result = self.supabase.table("golden_set_library").select(
                    "id, name, target_1, target_2, category, orr_pct, os_months, outcome_type, properties"
                ).neq("status", "rejected").or_(
                    f"description.ilike.%{term}%,category.ilike.%{term}%"
                ).limit(10).execute()

                for row in (result.data or []):
                    all_results.append({
                        "id": row["id"],
                        "name": row["name"],
                        "target_protein": row.get("target_1", "Unknown"),
                        "category": row.get("category"),
                        "clinical_score": self._score_golden_set_entry(row),
                        "combined_score": 0.8,
                        "similarity": 0.8
                    })

            # Deduplicate by target
            seen_targets = set()
            unique = []
            for ab in all_results:
                target = ab["target_protein"]
                if target not in seen_targets:
                    seen_targets.add(target)
                    unique.append(ab)

            return unique[:5]
        except Exception as e:
            logger.error(f"[navigator:step1] Golden set fallback error: {e}")
            return []

    async def _gemini_target_suggestion(self, disease_name: str) -> List[Dict]:
        """Gemini에게 질환에 맞는 타겟/항체 추천 요청"""
        try:
            model = get_gemini_model()
            prompt = f"""ADC(항체-약물접합체) 전문가로서 '{disease_name}' 치료를 위한
타겟 단백질과 항체 후보를 추천해주세요.

JSON 배열로 응답하세요:
[
  {{"name": "항체명", "target_protein": "타겟", "clinical_score": 0.0~1.0, "rationale": "추천 이유"}}
]

FDA 승인된 ADC를 우선 참고하되, 데이터가 없으면 '레퍼런스 부족'이라고 명시하세요."""

            response = await model.generate_content_async(prompt)
            content = response.text

            if "```json" in content:
                content = content.split("```json")[1].split("```")[0]
            elif "```" in content:
                content = content.split("```")[1].split("```")[0]

            suggestions = json.loads(content)
            results = []
            for s in suggestions[:5]:
                results.append({
                    "id": str(uuid.uuid4()),
                    "name": s.get("name", "AI-Suggested"),
                    "target_protein": s.get("target_protein", "Unknown"),
                    "clinical_score": float(s.get("clinical_score", 0.5)),
                    "combined_score": float(s.get("clinical_score", 0.5)),
                    "similarity": 0.6,
                    "source": "gemini_suggestion"
                })
            return results
        except Exception as e:
            logger.error(f"[navigator:step1] Gemini target suggestion error: {e}")
            return []

    async def _check_failure_cases(self, session_id: str, target_protein: str):
        """실패 사례 DB 체크 → 경고 발생"""
        try:
            result = self.supabase.table("golden_set_library").select(
                "name, target_1, outcome_type, failure_reason"
            ).neq("status", "rejected").eq(
                "outcome_type", "Failure"
            ).or_(
                f"target_1.ilike.%{target_protein}%,target_2.ilike.%{target_protein}%"
            ).execute()

            for failure in (result.data or []):
                logger.warning(
                    f"[navigator:step1] ⚠️ FAILURE CASE: {failure['name']} "
                    f"(target: {failure.get('target_1')}) - {failure.get('failure_reason', 'unknown')}"
                )
                await websocket_hub.stream_agent_log(
                    session_id, "warning",
                    f"⚠️ 실패 사례 감지: {failure['name']} - {failure.get('failure_reason', '')}",
                    emoji="⚠️", agent_name="librarian"
                )
        except Exception as e:
            logger.debug(f"[navigator:step1] Failure case check error: {e}")

    # =====================================================================
    # Step 2: Assembly via Alchemist Agent
    # =====================================================================

    async def _step2_assembly(
        self,
        session_id: str,
        target_protein: str,
        disease_name: str,
        antibody_candidates: List[AntibodyCandidate]
    ) -> tuple:
        """
        Step 2: AlchemistAgent.execute()를 호출하여
        실제 Gemini 기반 SMILES 생성 + Golden Set 참조로 조합을 설계합니다.
        """
        logger.info(f"[navigator:step2] Calling Alchemist for target: {target_protein}")

        # DesignSessionState 구성 (Alchemist가 읽는 필드)
        design_state: DesignSessionState = create_initial_state(
            session_id=session_id,
            user_id="navigator",
            session_type="denovo",
            tier="premium",
            target_antigen=target_protein,
            target_indication=disease_name,
            requested_dar=4,
            linker_preference="any",
            design_goal=f"Design optimal ADC for {disease_name} targeting {target_protein}"
        )

        # Alchemist 실제 호출
        alchemist_output: AgentOutput = await self.alchemist.execute(design_state)

        if alchemist_output.success and alchemist_output.data:
            primary_smiles = alchemist_output.data.get("primary_smiles", "")
            candidates = alchemist_output.data.get("candidates", [])
            golden_set_refs = alchemist_output.data.get("golden_set_refs", [])

            logger.info(
                f"[navigator:step2] Alchemist success: {len(candidates)} candidates, "
                f"SMILES length={len(primary_smiles)}, "
                f"golden_set_refs={len(golden_set_refs)}"
            )

            # Golden combination 구성
            golden_combo = self._build_golden_combination(
                antibody_candidates[0],
                primary_smiles,
                alchemist_output.data,
                target_protein
            )

            combined_smiles = primary_smiles
            return golden_combo, combined_smiles

        else:
            # Alchemist 실패 시 golden_set에서 직접 조합 조회
            logger.warning(f"[navigator:step2] Alchemist failed: {alchemist_output.error}, using golden_set fallback")
            return await self._golden_set_combination_fallback(
                target_protein, disease_name, antibody_candidates[0]
            )

    def _build_golden_combination(
        self,
        best_antibody: AntibodyCandidate,
        smiles: str,
        alchemist_data: Dict,
        target_protein: str
    ) -> GoldenCombination:
        """Alchemist 결과에서 GoldenCombination 구성"""
        # SMILES 분할 (linker.payload)
        smiles_parts = smiles.split(".") if smiles else ["", ""]
        linker_smiles = smiles_parts[0] if len(smiles_parts) > 1 else ""
        payload_smiles = smiles_parts[-1] if smiles_parts else ""

        return GoldenCombination(
            antibody=best_antibody,
            linker=LinkerSpec(
                id=str(uuid.uuid4()),
                smiles=linker_smiles,
                type=alchemist_data.get("linker_type", "cleavable"),
                cleavable=True
            ),
            payload=PayloadSpec(
                id=str(uuid.uuid4()),
                smiles=payload_smiles,
                class_name=alchemist_data.get("payload_class", "Unknown"),
                mechanism=alchemist_data.get("mechanism", None)
            ),
            dar=4,
            historical_performance=alchemist_data.get("historical_performance"),
            confidence_score=alchemist_data.get("confidence", 0.7)
        )

    async def _golden_set_combination_fallback(
        self,
        target_protein: str,
        disease_name: str,
        best_antibody: AntibodyCandidate
    ) -> tuple:
        """Golden set에서 직접 FDA 승인 ADC 조합 조회"""
        try:
            result = self.supabase.table("golden_set_library").select(
                "name, target_1, linker_smiles, payload_smiles, linker_type, dar, "
                "orr_pct, pfs_months, os_months, properties"
            ).neq("status", "rejected").or_(
                f"target_1.ilike.%{target_protein}%,target_2.ilike.%{target_protein}%"
            ).not_.is_("payload_smiles", "null").limit(5).execute()

            if result.data:
                best = result.data[0]
                props = best.get("properties") or {}

                linker_smiles = best.get("linker_smiles", "")
                payload_smiles = best.get("payload_smiles", "")
                combined = f"{linker_smiles}.{payload_smiles}" if linker_smiles and payload_smiles else (linker_smiles or payload_smiles)

                golden_combo = GoldenCombination(
                    antibody=best_antibody,
                    linker=LinkerSpec(
                        id=str(uuid.uuid4()),
                        smiles=linker_smiles,
                        type=best.get("linker_type", "cleavable"),
                        cleavable="cleavable" in (best.get("linker_type", "") or "").lower()
                    ),
                    payload=PayloadSpec(
                        id=str(uuid.uuid4()),
                        smiles=payload_smiles,
                        class_name=props.get("payload_class", "Unknown"),
                        mechanism=props.get("mechanism_of_action")
                    ),
                    dar=best.get("dar") or 4,
                    historical_performance={
                        "orr_pct": best.get("orr_pct"),
                        "pfs_months": best.get("pfs_months"),
                        "os_months": best.get("os_months"),
                        "source": "golden_set_library"
                    },
                    confidence_score=0.9
                )
                return golden_combo, combined

            # 전역 최고 성능 조합
            fallback = self.supabase.table("golden_set_library").select(
                "name, target_1, linker_smiles, payload_smiles, linker_type, dar, "
                "orr_pct, properties"
            ).neq("status", "rejected").not_.is_("payload_smiles", "null").order(
                "orr_pct", desc=True
            ).limit(1).execute()

            if fallback.data:
                fb = fallback.data[0]
                props = fb.get("properties") or {}
                ls = fb.get("linker_smiles", "")
                ps = fb.get("payload_smiles", "")
                combined = f"{ls}.{ps}" if ls and ps else (ls or ps)

                golden_combo = GoldenCombination(
                    antibody=best_antibody,
                    linker=LinkerSpec(id=str(uuid.uuid4()), smiles=ls,
                                     type=fb.get("linker_type", "cleavable")),
                    payload=PayloadSpec(id=str(uuid.uuid4()), smiles=ps,
                                       class_name=props.get("payload_class", "Unknown")),
                    dar=fb.get("dar") or 4,
                    historical_performance={"orr_pct": fb.get("orr_pct"), "source": "golden_set_global_best"},
                    confidence_score=0.6
                )
                return golden_combo, combined

            # 최종 fallback - Gemini에게 요청
            return await self._gemini_combination_fallback(target_protein, disease_name, best_antibody)

        except Exception as e:
            logger.error(f"[navigator:step2] Golden set fallback error: {e}")
            return await self._gemini_combination_fallback(target_protein, disease_name, best_antibody)

    async def _gemini_combination_fallback(
        self, target_protein: str, disease_name: str, best_antibody: AntibodyCandidate
    ) -> tuple:
        """Gemini에게 조합 추천 요청 (모든 DB 검색 실패 시)"""
        try:
            model = get_gemini_model()
            prompt = f"""ADC 전문가로서 {disease_name} 치료를 위한 {target_protein} 타겟 ADC 조합을 설계하세요.

JSON으로 응답:
{{
  "linker_smiles": "유효한 SMILES",
  "linker_type": "cleavable 또는 non-cleavable",
  "payload_smiles": "유효한 SMILES",
  "payload_class": "페이로드 분류명",
  "mechanism": "작용기전",
  "dar": 정수,
  "confidence": 0.0~1.0,
  "rationale": "설계 근거. 레퍼런스가 부족하면 명시할 것"
}}"""

            response = await model.generate_content_async(prompt)
            content = response.text
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0]
            elif "```" in content:
                content = content.split("```")[1].split("```")[0]

            data = json.loads(content)
            ls = data.get("linker_smiles", "")
            ps = data.get("payload_smiles", "")
            combined = f"{ls}.{ps}" if ls and ps else (ls or ps)

            golden_combo = GoldenCombination(
                antibody=best_antibody,
                linker=LinkerSpec(id=str(uuid.uuid4()), smiles=ls,
                                  type=data.get("linker_type", "cleavable")),
                payload=PayloadSpec(id=str(uuid.uuid4()), smiles=ps,
                                    class_name=data.get("payload_class", "AI-Designed"),
                                    mechanism=data.get("mechanism")),
                dar=data.get("dar", 4),
                historical_performance={"source": "gemini_design", "rationale": data.get("rationale", "")},
                confidence_score=data.get("confidence", 0.5)
            )
            return golden_combo, combined

        except Exception as e:
            logger.error(f"[navigator:step2] Gemini fallback error: {e}")
            # 진짜 최후의 수단 - 데이터 없음을 명시하는 빈 조합
            golden_combo = GoldenCombination(
                antibody=best_antibody,
                linker=LinkerSpec(id=str(uuid.uuid4()), smiles="", type="unknown"),
                payload=PayloadSpec(id=str(uuid.uuid4()), smiles="", class_name="Insufficient Data"),
                dar=4,
                historical_performance={"source": "no_data", "warning": "충분한 레퍼런스가 부족합니다"},
                confidence_score=0.1
            )
            return golden_combo, ""

    # =====================================================================
    # Step 3: Property Calculation via Coder Agent
    # =====================================================================

    async def _step3_property_calculation(
        self, session_id: str, smiles: str, target_protein: str, disease_name: str
    ) -> Dict[str, Any]:
        """
        Step 3: CoderAgent.execute()를 호출하여
        Sandbox(Docker/subprocess)에서 실제 RDKit 물성 계산을 수행합니다.
        """
        if not smiles or not smiles.strip():
            logger.warning("[navigator:step3] No SMILES available, skipping property calculation")
            return {"warning": "No SMILES available for calculation"}

        logger.info(f"[navigator:step3] Calling Coder for SMILES: {smiles[:50]}...")

        # DesignSessionState 구성 (Coder가 읽는 필드)
        design_state: DesignSessionState = create_initial_state(
            session_id=session_id,
            user_id="navigator",
            session_type="denovo",
            tier="premium",
            target_antigen=target_protein,
            target_indication=disease_name,
        )
        design_state["current_smiles"] = smiles

        # Coder 실제 호출
        coder_output: AgentOutput = await self.coder.execute(design_state)

        if coder_output.success:
            metrics = coder_output.data.get("metrics", {})
            logger.info(f"[navigator:step3] Coder success: {metrics}")
            return dict(metrics) if metrics else {}

        # Coder 실패 → Healer 호출 (최대 3회)
        logger.warning(f"[navigator:step3] Coder failed: {coder_output.error}, calling Healer")

        for attempt in range(3):
            design_state["healing_attempts"] = attempt + 1
            design_state["requires_healing"] = True
            design_state["last_error"] = coder_output.error
            design_state["last_code"] = coder_output.data.get("code", "")

            healer_output = await self.healer.execute(design_state)

            if healer_output.success:
                # Healer가 고친 코드로 Coder 재실행
                design_state["requires_healing"] = False
                coder_retry = await self.coder.execute(design_state)

                if coder_retry.success:
                    metrics = coder_retry.data.get("metrics", {})
                    logger.info(f"[navigator:step3] Coder success after healing (attempt {attempt + 1})")
                    return dict(metrics) if metrics else {}

                coder_output = coder_retry
            else:
                logger.warning(f"[navigator:step3] Healer failed at attempt {attempt + 1}")

        # Coder + Healer 모두 실패 → RDKit 직접 fallback
        logger.warning("[navigator:step3] All Coder+Healer attempts failed, using direct RDKit fallback")
        return self._direct_rdkit_calculation(smiles)

    def _direct_rdkit_calculation(self, smiles: str) -> Dict[str, Any]:
        """RDKit 직접 계산 (Coder+Healer 모두 실패 시 최후의 수단)"""
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen

            metrics = {}
            for smi in smiles.split("."):
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    mw = Descriptors.MolWt(mol)
                    logp = Crippen.MolLogP(mol)
                    metrics = {
                        "mw": round(mw, 2),
                        "logp": round(logp, 2),
                        "hbd": rdMolDescriptors.CalcNumHBD(mol),
                        "hba": rdMolDescriptors.CalcNumHBA(mol),
                        "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
                        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
                        "source": "direct_rdkit_fallback"
                    }
                    break  # 첫 번째 유효한 분자만

            return metrics
        except Exception as e:
            logger.error(f"[navigator:step3] Direct RDKit fallback error: {e}")
            return {"error": str(e), "source": "calculation_failed"}

    # =====================================================================
    # Step 4: Validation via Auditor Agent
    # =====================================================================

    async def _step4_validation(
        self,
        session_id: str,
        smiles: str,
        target_protein: str,
        disease_name: str,
        calculated_metrics: Dict[str, Any]
    ) -> tuple:
        """
        Step 4: AuditorAgent.execute()를 호출하여
        실제 PAINS/Lipinski/SA Score/독성 검증 + 실패 사례 체크를 수행합니다.
        """
        failure_warnings = []

        if not smiles or not smiles.strip():
            logger.warning("[navigator:step4] No SMILES, skipping validation")
            return False, [{"warning": "No SMILES available"}], []

        logger.info(f"[navigator:step4] Calling Auditor for validation")

        # DesignSessionState 구성 (Auditor가 읽는 필드)
        design_state: DesignSessionState = create_initial_state(
            session_id=session_id,
            user_id="navigator",
            session_type="denovo",
            tier="premium",
            target_antigen=target_protein,
            target_indication=disease_name,
        )
        design_state["current_smiles"] = smiles
        design_state["calculated_metrics"] = calculated_metrics

        # Auditor 실제 호출
        auditor_output: AgentOutput = await self.auditor.execute(design_state)

        validations = []
        physics_verified = False

        if auditor_output.success:
            decision = auditor_output.data.get("decision", {})
            chemistry = auditor_output.data.get("chemistry_validation", {})
            constraint = auditor_output.data.get("constraint_check", {})

            physics_verified = decision.get("approved", False)

            validations.append({
                "source": "auditor_agent",
                "decision": decision.get("action", "unknown"),
                "approved": decision.get("approved", False),
                "confidence": decision.get("confidence", 0),
                "reasoning": decision.get("reasoning", ""),
                "suggestions": decision.get("suggestions", []),
                "key_concerns": decision.get("key_concerns", []),
                "chemistry": chemistry,
                "constraints": constraint,
                "risk_score": auditor_output.data.get("risk_score", 0),
            })
        else:
            logger.warning(f"[navigator:step4] Auditor returned failure: {auditor_output.error}")
            validations.append({
                "source": "auditor_agent",
                "error": auditor_output.error,
                "approved": False,
            })

        # 실패 사례 교차 체크 (Auditor 결과와 별도로)
        failure_warnings = await self._check_failure_similarity(session_id, smiles, target_protein)

        return physics_verified, validations, failure_warnings

    async def _check_failure_similarity(
        self, session_id: str, smiles: str, target_protein: str
    ) -> List[str]:
        """실패 사례와의 유사성 체크"""
        warnings = []
        try:
            failures = self.supabase.table("golden_set_library").select(
                "name, target_1, failure_reason, properties"
            ).neq("status", "rejected").eq(
                "outcome_type", "Failure"
            ).execute()

            for failure in (failures.data or []):
                props = failure.get("properties") or {}
                failure_name = failure.get("name", "")

                # 타겟 일치 체크
                if target_protein and failure.get("target_1", "").lower() == target_protein.lower():
                    warning = f"동일 타겟 실패 사례: {failure_name} - {failure.get('failure_reason', '')}"
                    warnings.append(warning)
                    await websocket_hub.stream_agent_log(
                        session_id, "warning", f"⚠️ {warning}",
                        emoji="⚠️", agent_name="auditor"
                    )

                # PBD dimer 간독성 체크
                payload_class = props.get("payload_class", "")
                if "PBD" in payload_class.upper() and "PBD" in smiles.upper():
                    warning = f"PBD dimer 간독성 위험 (참고: {failure_name})"
                    warnings.append(warning)
                    await websocket_hub.stream_agent_log(
                        session_id, "warning", f"🚨 {warning}",
                        emoji="🚨", agent_name="auditor"
                    )

        except Exception as e:
            logger.debug(f"[navigator:step4] Failure similarity check error: {e}")

        return warnings

    # =====================================================================
    # Step 5: Virtual Clinical Trial (Gemini + Golden Set Reference)
    # =====================================================================

    async def _step5_virtual_trial(
        self,
        session_id: str,
        golden_combo: GoldenCombination,
        metrics: Dict[str, Any],
        disease_name: str,
        target_protein: str,
    ) -> VirtualTrialResult:
        """
        Step 5: 실제 golden_set 임상 데이터 기반 + Gemini AI 예측으로
        가상 임상 결과를 생성합니다.
        """
        logger.info(f"[navigator:step5] Virtual trial for {target_protein}/{disease_name}")

        # 1. Golden set에서 동일 타겟의 실제 임상 데이터 조회
        reference_data = await self._get_clinical_references(target_protein, disease_name)

        # 2. Gemini 기반 예측
        gemini_prediction = await self._gemini_clinical_prediction(
            target_protein, disease_name, golden_combo, metrics, reference_data
        )

        # 3. 레퍼런스 데이터와 Gemini 예측 결합
        predicted_orr = gemini_prediction.get("predicted_orr", 0)
        predicted_pfs = gemini_prediction.get("predicted_pfs_months", 0)
        predicted_os = gemini_prediction.get("predicted_os_months", 0)
        confidence = gemini_prediction.get("confidence", 0.5)

        # 4. PK/PD 시뮬레이션 (2-구획 모델 + 종양 성장)
        pk_data = self._simulate_pk(golden_combo)
        tumor_data = self._simulate_tumor_growth(predicted_orr)

        result = VirtualTrialResult(
            predicted_orr=predicted_orr,
            predicted_pfs_months=predicted_pfs,
            predicted_os_months=predicted_os,
            pk_data=pk_data,
            tumor_data=tumor_data,
            patient_population=disease_name,
            confidence=confidence
        )

        await websocket_hub.stream_agent_log(
            session_id, "success",
            f"임상 예측 완료: ORR={predicted_orr:.1f}%, PFS={predicted_pfs:.1f}mo, "
            f"신뢰도={confidence:.0%}",
            emoji="📊", agent_name="clinical"
        )

        return result

    async def _get_clinical_references(self, target_protein: str, disease_name: str) -> List[Dict]:
        """Golden set에서 실제 임상 데이터 조회"""
        try:
            search_terms = self._get_search_terms(disease_name)
            all_refs = []

            # 타겟 기반 검색
            result = self.supabase.table("golden_set_library").select(
                "name, target_1, orr_pct, pfs_months, os_months, outcome_type, properties"
            ).neq("status", "rejected").or_(
                f"target_1.ilike.%{target_protein}%,target_2.ilike.%{target_protein}%"
            ).execute()

            for row in (result.data or []):
                all_refs.append(row)

            # 질환 기반 추가 검색
            for term in search_terms[:2]:
                result2 = self.supabase.table("golden_set_library").select(
                    "name, target_1, orr_pct, pfs_months, os_months, outcome_type, properties"
                ).neq("status", "rejected").or_(
                    f"description.ilike.%{term}%,category.ilike.%{term}%"
                ).execute()

                for row in (result2.data or []):
                    if row not in all_refs:
                        all_refs.append(row)

            logger.info(f"[navigator:step5] Found {len(all_refs)} clinical references")
            return all_refs

        except Exception as e:
            logger.error(f"[navigator:step5] Clinical reference error: {e}")
            return []

    async def _gemini_clinical_prediction(
        self,
        target_protein: str,
        disease_name: str,
        golden_combo: GoldenCombination,
        metrics: Dict[str, Any],
        reference_data: List[Dict]
    ) -> Dict[str, Any]:
        """Gemini에게 실제 레퍼런스 기반 임상 예측 요청"""
        try:
            model = get_gemini_model()

            # 레퍼런스 요약
            ref_summary = ""
            for ref in reference_data[:8]:
                ref_summary += (
                    f"- {ref.get('name', 'Unknown')}: "
                    f"Target={ref.get('target_1', 'N/A')}, "
                    f"ORR={ref.get('orr_pct', 'N/A')}%, "
                    f"PFS={ref.get('pfs_months', 'N/A')}mo, "
                    f"OS={ref.get('os_months', 'N/A')}mo, "
                    f"Outcome={ref.get('outcome_type', 'N/A')}\n"
                )

            if not ref_summary:
                ref_summary = "레퍼런스 데이터 없음. 일반적인 ADC 임상 결과 기반으로 보수적으로 예측하세요."

            # 물성 요약
            metrics_summary = json.dumps(metrics, indent=2) if metrics else "물성 데이터 없음"

            prompt = f"""ADC 임상 전문가로서 아래 설계된 ADC의 가상 임상 결과를 예측해주세요.

**설계 ADC 정보**
- 타겟: {target_protein}
- 적응증: {disease_name}
- 항체: {golden_combo.antibody.name}
- 링커: {golden_combo.linker.type}
- 페이로드: {golden_combo.payload.class_name}
- DAR: {golden_combo.dar}

**물성 정보**
{metrics_summary}

**레퍼런스 임상 데이터 (FDA 승인 ADC)**
{ref_summary}

레퍼런스가 부족하면 '충분한 레퍼런스가 부족함'을 명시하고 보수적으로 예측하세요.
실패 사례가 있으면 반드시 반영하세요.

JSON으로 응답:
{{
  "predicted_orr": 0~100 (숫자),
  "predicted_pfs_months": 양수 (숫자),
  "predicted_os_months": 양수 (숫자),
  "confidence": 0.0~1.0,
  "reasoning": "예측 근거 (어떤 레퍼런스를 참고했는지 명시)",
  "limitations": "예측 제한사항"
}}"""

            response = await model.generate_content_async(prompt)
            content = response.text

            if "```json" in content:
                content = content.split("```json")[1].split("```")[0]
            elif "```" in content:
                content = content.split("```")[1].split("```")[0]

            prediction = json.loads(content)

            # 값 검증
            predicted_orr = max(0, min(100, float(prediction.get("predicted_orr", 0))))
            predicted_pfs = max(0, float(prediction.get("predicted_pfs_months", 0)))
            predicted_os = max(0, float(prediction.get("predicted_os_months", 0)))
            confidence = max(0, min(1.0, float(prediction.get("confidence", 0.5))))

            logger.info(
                f"[navigator:step5] Gemini prediction: "
                f"ORR={predicted_orr}%, PFS={predicted_pfs}mo, "
                f"confidence={confidence}"
            )

            return {
                "predicted_orr": predicted_orr,
                "predicted_pfs_months": predicted_pfs,
                "predicted_os_months": predicted_os,
                "confidence": confidence,
                "reasoning": prediction.get("reasoning", ""),
                "limitations": prediction.get("limitations", "")
            }

        except Exception as e:
            logger.error(f"[navigator:step5] Gemini prediction error: {e}")

            # 레퍼런스 기반 통계적 추정 (Gemini 실패 시)
            if reference_data:
                orrs = [r.get("orr_pct", 0) for r in reference_data if r.get("orr_pct")]
                pfss = [r.get("pfs_months", 0) for r in reference_data if r.get("pfs_months")]
                oss = [r.get("os_months", 0) for r in reference_data if r.get("os_months")]

                return {
                    "predicted_orr": sum(orrs) / len(orrs) if orrs else 0,
                    "predicted_pfs_months": sum(pfss) / len(pfss) if pfss else 0,
                    "predicted_os_months": sum(oss) / len(oss) if oss else 0,
                    "confidence": 0.4,
                    "reasoning": f"통계적 추정 (레퍼런스 {len(reference_data)}건 기반)",
                    "limitations": "Gemini AI 예측 실패, 단순 통계 기반"
                }

            return {
                "predicted_orr": 0,
                "predicted_pfs_months": 0,
                "predicted_os_months": 0,
                "confidence": 0.1,
                "reasoning": "충분한 레퍼런스가 부족합니다",
                "limitations": "데이터 및 AI 예측 모두 실패"
            }

    def _simulate_pk(self, golden_combo: GoldenCombination) -> List[Dict[str, float]]:
        """2-구획 PK 시뮬레이션 (DAR/payload class 기반 파라미터)"""
        # Payload class별 PK 파라미터 (문헌 기반)
        PK_PARAMS = {
            "MMAE": {"clearance": 0.04, "vd": 6.0, "half_life": 72},
            "DXd": {"clearance": 0.03, "vd": 5.5, "half_life": 144},
            "DM1": {"clearance": 0.03, "vd": 5.0, "half_life": 96},
            "SN-38": {"clearance": 0.05, "vd": 7.0, "half_life": 120},
            "PBD": {"clearance": 0.02, "vd": 4.0, "half_life": 48},
            "MMAF": {"clearance": 0.035, "vd": 5.5, "half_life": 80},
        }

        payload_class = (golden_combo.payload.class_name or "").upper()
        pk = None
        for key, params in PK_PARAMS.items():
            if key.upper() in payload_class:
                pk = params
                break
        if not pk:
            pk = {"clearance": 0.04, "vd": 6.0, "half_life": 96}

        # DAR 보정
        dar = golden_combo.dar or 4
        cl = pk["clearance"] * (1 + (dar - 4) * 0.05)
        vd = pk["vd"]
        half_life = pk["half_life"] * (1 - (dar - 4) * 0.03)

        # 2-구획 모델 파라미터
        dose = 3.6  # mg/kg
        alpha = 0.693 / (half_life * 0.3)
        beta = 0.693 / half_life
        cmax = dose * 1000 / vd

        pk_data = []
        for t in range(0, 505, 12):
            conc = cmax * (0.6 * math.exp(-alpha * t) + 0.4 * math.exp(-beta * t))
            free_payload = conc * cl * 0.01 * (1 - math.exp(-0.01 * t))
            pk_data.append({
                "time_hours": t,
                "concentration": round(max(0, conc), 3),
                "free_payload": round(max(0, free_payload), 3)
            })

        return pk_data

    def _simulate_tumor_growth(self, predicted_orr: float) -> List[Dict[str, float]]:
        """종양 성장 억제 시뮬레이션"""
        initial_volume = 200
        doubling_time = 15  # days
        effect = min(0.95, (predicted_orr or 0) / 100 * 0.9)

        tumor_data = []
        for day in range(0, 43, 3):
            control = initial_volume * (2 ** (day / doubling_time))
            treated = initial_volume * (2 ** (day / doubling_time)) * (1 - effect * (1 - math.exp(-0.1 * day)))
            tumor_data.append({
                "day": day,
                "treated": round(max(initial_volume * 0.3, treated), 1),
                "control": round(control, 1)
            })

        return tumor_data

    # =====================================================================
    # Utility Methods
    # =====================================================================

    def _get_search_terms(self, disease_name: str) -> List[str]:
        """질환명에서 검색어 목록 생성"""
        terms = [disease_name]

        # 주요 키워드 추출
        keywords = disease_name.lower().split()
        for keyword in keywords:
            if len(keyword) > 3 and keyword not in ["cancer", "disease", "syndrome"]:
                terms.append(keyword)

        # Subtype expansion
        subtype_map = {
            "breast cancer": ["breast", "HER2", "TNBC", "HR+"],
            "lung cancer": ["lung", "NSCLC", "SCLC"],
            "gastric cancer": ["gastric", "stomach", "GEJ"],
            "colorectal cancer": ["colorectal", "colon", "CRC"],
            "pancreatic cancer": ["pancreatic", "pancreas"],
            "bladder cancer": ["bladder", "urothelial"],
            "ovarian cancer": ["ovarian", "ovary"],
        }

        for key, expansions in subtype_map.items():
            if key in disease_name.lower():
                terms.extend(expansions)
                break

        return list(dict.fromkeys(terms))  # deduplicate preserving order

    def _score_golden_set_entry(self, row: Dict) -> float:
        """Golden set 항목의 임상 점수 계산"""
        score = 0.0
        orr = row.get("orr_pct", 0) or 0
        score += (orr / 100) * 0.35

        pfs = row.get("pfs_months", 0) or 0
        score += min(pfs / 24, 1.0) * 0.25

        os_m = row.get("os_months", 0) or 0
        score += min(os_m / 36, 1.0) * 0.20

        outcome = row.get("outcome_type", "")
        if outcome in ("Approved", "Success"):
            score += 0.15
        elif outcome == "Failure":
            score -= 0.1

        return round(max(0, score), 4)

    # =====================================================================
    # DB Operations
    # =====================================================================

    async def _create_db_session(self, state: NavigatorState):
        try:
            self.supabase.table("navigator_sessions").upsert({
                "id": state.session_id,
                "user_id": state.user_id,
                "disease_name": state.disease_name,
                "status": "running",
                "current_step": 0,
                "total_steps": state.total_steps,
                "created_at": datetime.utcnow().isoformat()
            }).execute()
        except Exception as e:
            logger.error(f"[navigator] DB session create error: {e}")

    async def _update_db_step(self, state: NavigatorState, step: int, data: Dict = None):
        try:
            update = {
                "current_step": step,
                "updated_at": datetime.utcnow().isoformat()
            }
            if data:
                for key, value in data.items():
                    update[key] = value
            self.supabase.table("navigator_sessions").update(update).eq(
                "id", state.session_id
            ).execute()
        except Exception as e:
            logger.error(f"[navigator] DB step update error: {e}")

    async def _update_db_error(self, state: NavigatorState, error: str):
        try:
            self.supabase.table("navigator_sessions").update({
                "status": "failed",
                "error_message": error,
                "updated_at": datetime.utcnow().isoformat()
            }).eq("id", state.session_id).execute()
        except Exception as e:
            logger.error(f"[navigator] DB error update error: {e}")

    async def _broadcast_step(self, state: NavigatorState, step: int, message: str):
        try:
            await websocket_hub.broadcast_agent_status(
                session_id=state.session_id,
                agent="navigator",
                status="running",
                message=message,
                step=step
            )
            logger.info(f"[navigator:step{step}] {message}")
        except Exception as e:
            logger.debug(f"[navigator] Broadcast error: {e}")


# ============================================================================
# Public API
# ============================================================================

async def run_one_click_navigator(
    disease_name: str,
    session_id: Optional[str] = None,
    user_id: Optional[str] = None
) -> NavigatorResult:
    """Navigator 실행 진입점"""
    orchestrator = NavigatorOrchestrator()
    return await orchestrator.run_one_click_pipeline(
        disease_name=disease_name,
        session_id=session_id,
        user_id=user_id
    )
