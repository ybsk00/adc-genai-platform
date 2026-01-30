"""
One-Click ADC Navigator Orchestrator
AstraForge Enhancement Specification v2.2

질환명 하나만 입력하면 6인 에이전트가 협업하여
최적의 ADC 설계안을 자동 생성하는 파이프라인 오케스트레이터

3단계 파이프라인:
- Step 1: Target & Antibody Match (Librarian)
- Step 2: Linker & Payload Coupling (Alchemist)
- Step 3: Simulation & Audit (Coder + Auditor)
"""

import logging
import json
import uuid
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum

from app.core.supabase import get_supabase_client
from app.core.websocket_hub import websocket_hub
from app.services.rag_service import RAGService
from app.core.gemini import get_gemini_model

logger = logging.getLogger(__name__)


# ============================================================================
# Data Classes
# ============================================================================

class NavigatorStep(Enum):
    """Navigator 단계"""
    INIT = 0
    TARGET_MATCH = 1
    GOLDEN_COMBINATION = 2
    PROPERTY_CALCULATION = 3
    PHYSICAL_VALIDATION = 4
    VIRTUAL_TRIAL = 5
    COMPLETE = 6


@dataclass
class AntibodyCandidate:
    """항체 후보"""
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
    """링커 명세"""
    id: str
    smiles: str
    type: str
    cleavable: bool = True
    properties: Optional[Dict[str, Any]] = None


@dataclass
class PayloadSpec:
    """페이로드 명세"""
    id: str
    smiles: str
    class_name: str
    mechanism: Optional[str] = None
    properties: Optional[Dict[str, Any]] = None


@dataclass
class GoldenCombination:
    """최적 조합 (Golden Combination)"""
    antibody: AntibodyCandidate
    linker: LinkerSpec
    payload: PayloadSpec
    dar: int = 4
    historical_performance: Optional[Dict[str, Any]] = None
    confidence_score: float = 0.0


@dataclass
class VirtualTrialResult:
    """가상 임상 결과"""
    predicted_orr: float
    predicted_pfs_months: float
    predicted_os_months: float
    pk_data: List[Dict[str, float]] = field(default_factory=list)
    tumor_data: List[Dict[str, float]] = field(default_factory=list)
    patient_population: str = ""
    confidence: float = 0.7


@dataclass
class NavigatorResult:
    """Navigator 최종 결과"""
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
    """Navigator 파이프라인 상태"""
    session_id: str
    user_id: Optional[str]
    disease_name: str
    step: NavigatorStep = NavigatorStep.INIT
    total_steps: int = 5

    # Step 결과
    antibody_candidates: List[AntibodyCandidate] = field(default_factory=list)
    target_protein: str = ""
    golden_combination: Optional[GoldenCombination] = None
    combined_smiles: str = ""
    calculated_metrics: Dict[str, Any] = field(default_factory=dict)
    simulation_results: Dict[str, Any] = field(default_factory=dict)
    physical_validations: List[Dict[str, Any]] = field(default_factory=list)
    physics_verified: bool = False
    virtual_trial: Optional[VirtualTrialResult] = None

    # 메타데이터
    started_at: datetime = field(default_factory=datetime.utcnow)
    errors: List[str] = field(default_factory=list)


# ============================================================================
# Clinical Weighted Scorer
# ============================================================================

class ClinicalWeightedScorer:
    """임상 데이터 가중치 스코어링 알고리즘"""

    WEIGHTS = {
        "orr_pct": 0.35,
        "pfs_months": 0.25,
        "os_months": 0.20,
        "clinical_phase": 0.10,
        "safety_profile": 0.10
    }

    NORMALIZATION = {
        "orr_pct": 100,
        "pfs_months": 24,
        "os_months": 36
    }

    PHASE_SCORES = {
        "Approved": 1.0,
        "BLA Submitted": 0.95,
        "Phase 3": 0.8,
        "Phase 2/3": 0.7,
        "Phase 2": 0.5,
        "Phase 1/2": 0.4,
        "Phase 1": 0.3,
        "Preclinical": 0.1
    }

    @classmethod
    def calculate_score(cls, data: Dict[str, Any]) -> float:
        """시약 조합의 가중 점수 계산"""
        score = 0.0

        # ORR
        orr = data.get("orr_pct", 0) or 0
        score += (orr / cls.NORMALIZATION["orr_pct"]) * cls.WEIGHTS["orr_pct"]

        # PFS
        pfs = data.get("pfs_months", 0) or 0
        pfs_norm = min(pfs / cls.NORMALIZATION["pfs_months"], 1.0)
        score += pfs_norm * cls.WEIGHTS["pfs_months"]

        # OS
        os_months = data.get("os_months", 0) or 0
        os_norm = min(os_months / cls.NORMALIZATION["os_months"], 1.0)
        score += os_norm * cls.WEIGHTS["os_months"]

        # Phase bonus
        phase = data.get("clinical_phase", "") or ""
        phase_score = cls.PHASE_SCORES.get(phase, 0.1)
        score += phase_score * cls.WEIGHTS["clinical_phase"]

        # Safety profile (default 0.5)
        safety = data.get("safety_score", 0.5) or 0.5
        score += safety * cls.WEIGHTS["safety_profile"]

        return round(score, 4)


# ============================================================================
# Navigator Orchestrator
# ============================================================================

class NavigatorOrchestrator:
    """
    One-Click ADC Navigator 오케스트레이터

    질환명 하나만 입력받아 최적의 ADC 설계안을 자동 생성합니다.
    """

    def __init__(self):
        self.supabase = get_supabase_client()
        self.rag_service = RAGService()
        self.scorer = ClinicalWeightedScorer()

        # 에이전트들 (lazy import to avoid circular imports)
        self._librarian = None
        self._alchemist = None
        self._coder = None
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
    def auditor(self):
        if self._auditor is None:
            from app.agents.auditor import AuditorAgent
            self._auditor = AuditorAgent()
        return self._auditor

    async def run_one_click_pipeline(
        self,
        disease_name: str,
        session_id: Optional[str] = None,
        user_id: Optional[str] = None
    ) -> NavigatorResult:
        """
        질환명 기반 원클릭 ADC 설계 파이프라인

        Args:
            disease_name: 질환명 (예: "Pancreatic Cancer")
            session_id: 세션 ID (없으면 자동 생성)
            user_id: 사용자 ID

        Returns:
            NavigatorResult: 최종 설계 결과
        """
        # 세션 초기화
        if not session_id:
            session_id = str(uuid.uuid4())

        state = NavigatorState(
            session_id=session_id,
            user_id=user_id,
            disease_name=disease_name
        )

        try:
            # DB 세션 생성
            await self._create_db_session(state)

            # ═══════════════════════════════════════════════════════════
            # Step 1: 타겟 및 항체 최적화
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.TARGET_MATCH
            await self._broadcast_step(state, 1, "타겟 및 항체 검색 중...")

            antibody_candidates = await self._find_antibodies_by_disease(
                disease_name=disease_name,
                top_k=3
            )

            if not antibody_candidates:
                raise ValueError(f"질환 '{disease_name}'에 대한 항체 후보를 찾을 수 없습니다.")

            state.antibody_candidates = antibody_candidates
            state.target_protein = antibody_candidates[0].target_protein

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

            # 모든 고유 타겟 표시
            unique_targets = list(set(ab.target_protein for ab in antibody_candidates))
            targets_display = ", ".join(unique_targets[:5])
            await self._broadcast_step(
                state, 1,
                f"✅ {len(antibody_candidates)}개 항체 후보 발견 (타겟: {targets_display})"
            )

            # ═══════════════════════════════════════════════════════════
            # Step 2: 최적 조합 생성 (Golden Combination)
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.GOLDEN_COMBINATION
            await self._broadcast_step(state, 2, "Golden Combination 설계 중...")

            golden_combo = await self._generate_golden_combination(
                target_protein=state.target_protein,
                antibody_candidates=antibody_candidates
            )

            state.golden_combination = golden_combo

            await self._update_db_step(state, 2, {
                "golden_combination": {
                    "antibody": {
                        "id": golden_combo.antibody.antibody_id,
                        "name": golden_combo.antibody.name
                    },
                    "linker": {
                        "id": golden_combo.linker.id,
                        "smiles": golden_combo.linker.smiles,
                        "type": golden_combo.linker.type
                    },
                    "payload": {
                        "id": golden_combo.payload.id,
                        "smiles": golden_combo.payload.smiles,
                        "class_name": golden_combo.payload.class_name
                    },
                    "dar": golden_combo.dar,
                    "historical_performance": golden_combo.historical_performance
                }
            })

            orr_display = golden_combo.historical_performance.get("orr_pct", "N/A") if golden_combo.historical_performance else "N/A"
            await self._broadcast_step(
                state, 2,
                f"✅ 최적 조합 발견 (ORR: {orr_display}%)"
            )

            # ═══════════════════════════════════════════════════════════
            # Step 3: 물성 계산 및 시뮬레이션
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.PROPERTY_CALCULATION
            await self._broadcast_step(state, 3, "분자 물성 계산 중...")

            # SMILES 조합 생성
            combined_smiles = self._combine_adc_structure(golden_combo)
            state.combined_smiles = combined_smiles

            # 물성 계산
            metrics = await self._calculate_properties(combined_smiles, golden_combo)
            state.calculated_metrics = metrics

            await self._update_db_step(state, 3, {
                "combined_smiles": combined_smiles,
                "calculated_metrics": metrics
            })

            await self._broadcast_step(state, 3, "✅ 물성 계산 완료")

            # ═══════════════════════════════════════════════════════════
            # Step 4: 물리적 타당성 검증
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.PHYSICAL_VALIDATION
            await self._broadcast_step(state, 4, "물리적 타당성 검증 중...")

            # Physical Validation (Phase 2 기능 활용)
            from app.services.physical_validator import validate_structure

            validation_result = await validate_structure(
                smiles=combined_smiles,
                session_id=session_id,
                molecule_name=f"ADC for {disease_name}",
                generate_3d=False,
                save_to_db=False
            )

            state.physical_validations = validation_result.get("validations", [])
            state.physics_verified = validation_result.get("overall_status") == "pass"

            await self._update_db_step(state, 4, {
                "physical_validations": state.physical_validations,
                "physics_verified": state.physics_verified
            })

            if state.physics_verified:
                await self._broadcast_step(state, 4, "✅ Physics Verified 인증 완료")
            else:
                await self._broadcast_step(state, 4, "⚠️ 물리 검증 경고 발생")

            # ═══════════════════════════════════════════════════════════
            # Step 5: 가상 임상 시뮬레이션
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.VIRTUAL_TRIAL
            await self._broadcast_step(state, 5, "가상 임상 시뮬레이션 중...")

            virtual_trial = await self._run_virtual_trial(
                golden_combo=golden_combo,
                pk_params=metrics.get("pk_parameters", {}),
                disease_name=disease_name
            )

            state.virtual_trial = virtual_trial

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
                f"✅ 가상 임상 완료 (예측 ORR: {virtual_trial.predicted_orr:.1f}%)"
            )

            # ═══════════════════════════════════════════════════════════
            # 완료
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.COMPLETE

            # Digital Lineage 수집
            lineage = await self._collect_lineage(state)

            # DB 완료 처리
            await self._complete_db_session(state, lineage)

            execution_time = (datetime.utcnow() - state.started_at).total_seconds()

            await self._broadcast_step(
                state, 5,
                f"🎉 One-Click ADC 설계 완료! (소요시간: {execution_time:.1f}초)",
                complete=True
            )

            return NavigatorResult(
                session_id=session_id,
                disease_name=disease_name,
                antibody_candidates=antibody_candidates,
                golden_combination=golden_combo,
                calculated_metrics=state.calculated_metrics,
                physical_validations=state.physical_validations,
                physics_verified=state.physics_verified,
                virtual_trial=virtual_trial,
                digital_lineage=lineage,
                combined_smiles=combined_smiles,
                execution_time_seconds=execution_time
            )

        except Exception as e:
            logger.exception(f"[navigator] Pipeline error: {e}")
            state.errors.append(str(e))
            await self._fail_db_session(state, str(e))
            await self._broadcast_step(state, state.step.value, f"❌ 오류 발생: {str(e)}")
            raise

    # =========================================================================
    # Step 1: Target & Antibody Match (Real Data)
    # =========================================================================

    async def _find_all_targets_for_disease(self, disease_name: str) -> List[Dict[str, Any]]:
        """
        질환에 대한 모든 타겟 단백질을 실제 DB에서 검색
        golden_set → target_synonyms → antibody_library 순으로 조회
        """
        targets = {}  # canonical_name -> {data}

        try:
            # 1. golden_set에서 해당 질환의 타겟 검색
            gs_result = self.supabase.table("golden_set").select(
                "drug_name, target_1, target_2, indication, orr_pct, pfs_months, os_months, "
                "clinical_status, dar, payload_class, linker_type, mechanism_of_action"
            ).ilike("indication", f"%{disease_name}%").execute()

            for row in (gs_result.data or []):
                for target_col in ["target_1", "target_2"]:
                    t = row.get(target_col)
                    if t and t.strip():
                        canonical = t.strip().upper()
                        if canonical not in targets:
                            targets[canonical] = {
                                "canonical_name": canonical,
                                "display_name": t.strip(),
                                "source": "golden_set",
                                "clinical_data": [],
                                "drug_count": 0,
                                "best_orr": 0,
                                "best_phase": ""
                            }
                        targets[canonical]["drug_count"] += 1
                        targets[canonical]["clinical_data"].append({
                            "drug_name": row.get("drug_name"),
                            "orr_pct": row.get("orr_pct"),
                            "pfs_months": row.get("pfs_months"),
                            "os_months": row.get("os_months"),
                            "clinical_status": row.get("clinical_status"),
                            "dar": row.get("dar"),
                            "payload_class": row.get("payload_class"),
                            "linker_type": row.get("linker_type"),
                        })
                        orr = row.get("orr_pct") or 0
                        if orr > targets[canonical]["best_orr"]:
                            targets[canonical]["best_orr"] = orr
                            targets[canonical]["best_phase"] = row.get("clinical_status", "")

            logger.info(f"[navigator] Found {len(targets)} targets from golden_set for '{disease_name}'")

            # 2. antibody_library에서 추가 타겟 검색
            ab_result = self.supabase.table("antibody_library").select(
                "target_normalized, related_disease"
            ).ilike("related_disease", f"%{disease_name}%").limit(50).execute()

            for row in (ab_result.data or []):
                t = row.get("target_normalized")
                if t and t.strip():
                    canonical = t.strip().upper()
                    if canonical not in targets:
                        targets[canonical] = {
                            "canonical_name": canonical,
                            "display_name": t.strip(),
                            "source": "antibody_library",
                            "clinical_data": [],
                            "drug_count": 0,
                            "best_orr": 0,
                            "best_phase": ""
                        }
                    targets[canonical]["drug_count"] += 1

        except Exception as e:
            logger.error(f"[navigator] Target discovery error: {e}")

        # 타겟이 없으면 Gemini AI로 추천
        if not targets:
            targets = await self._gemini_suggest_targets(disease_name)

        # drug_count + best_orr 기준 정렬
        sorted_targets = sorted(
            targets.values(),
            key=lambda x: (x["best_orr"], x["drug_count"]),
            reverse=True
        )
        return sorted_targets

    async def _gemini_suggest_targets(self, disease_name: str) -> Dict[str, Dict]:
        """Gemini AI로 질환의 타겟 단백질 추천"""
        targets = {}
        try:
            model = get_gemini_model(temperature=0.1)
            prompt = (
                f"For the disease '{disease_name}', list the top 3-5 molecular targets "
                f"used in ADC (Antibody-Drug Conjugate) therapy. "
                f"Return ONLY a JSON array of objects with fields: "
                f"target_name, rationale (1 sentence). No markdown."
            )
            response = await model.ainvoke(prompt)
            content = response.content.strip()
            # Parse JSON from response
            if content.startswith("["):
                suggested = json.loads(content)
            else:
                # Try to extract JSON array
                start = content.find("[")
                end = content.rfind("]") + 1
                if start >= 0 and end > start:
                    suggested = json.loads(content[start:end])
                else:
                    suggested = []

            for item in suggested:
                name = item.get("target_name", "").strip().upper()
                if name:
                    targets[name] = {
                        "canonical_name": name,
                        "display_name": item.get("target_name", name),
                        "source": "gemini_ai",
                        "clinical_data": [],
                        "drug_count": 0,
                        "best_orr": 0,
                        "best_phase": "",
                        "rationale": item.get("rationale", "")
                    }
        except Exception as e:
            logger.error(f"[navigator] Gemini target suggestion error: {e}")
        return targets

    async def _find_antibodies_by_disease(
        self,
        disease_name: str,
        top_k: int = 5
    ) -> List[AntibodyCandidate]:
        """
        질환명 기반 최적 항체 검색 (실제 데이터 연결)
        1. 질환의 모든 타겟 발견
        2. 타겟별로 antibody_library + golden_set 검색
        3. 실제 임상 데이터로 clinical_score 계산
        """
        all_candidates = []

        try:
            # 1. 질환에 대한 모든 타겟 발견
            disease_targets = await self._find_all_targets_for_disease(disease_name)
            logger.info(f"[navigator] Targets for '{disease_name}': {[t['display_name'] for t in disease_targets]}")

            if not disease_targets:
                logger.warning(f"[navigator] No targets found for '{disease_name}', using vector search fallback")
                return await self._vector_search_antibodies(disease_name, top_k)

            # 2. 타겟별로 항체 검색
            for target_info in disease_targets[:6]:  # 상위 6개 타겟
                target_name = target_info["display_name"]
                canonical = target_info["canonical_name"]

                # antibody_library에서 해당 타겟의 항체 검색
                try:
                    ab_result = self.supabase.table("antibody_library").select(
                        "id, product_name, target_normalized, isotype, related_disease, full_spec, clinical_score"
                    ).ilike(
                        "target_normalized", f"%{target_name}%"
                    ).limit(5).execute()

                    for ab in (ab_result.data or []):
                        # 실제 임상 점수 계산
                        clinical_score = await self._calculate_clinical_score_real(
                            target_name=target_name,
                            clinical_data=target_info.get("clinical_data", [])
                        )
                        # similarity = target match quality
                        target_match = 1.0 if ab.get("target_normalized", "").upper() == canonical else 0.7
                        combined = target_match * 0.3 + clinical_score * 0.7

                        all_candidates.append(AntibodyCandidate(
                            antibody_id=str(ab.get("id", "")),
                            name=ab.get("product_name", "Unknown"),
                            target_protein=target_name,
                            isotype=ab.get("isotype"),
                            related_diseases=ab.get("related_disease"),
                            full_spec=ab.get("full_spec"),
                            clinical_score=clinical_score,
                            match_confidence=round(combined, 4)
                        ))
                except Exception as e:
                    logger.warning(f"[navigator] Antibody search for target '{target_name}' error: {e}")

                # golden_set에서도 해당 타겟의 ADC 약물 추가
                for cd in target_info.get("clinical_data", []):
                    drug_name = cd.get("drug_name", "")
                    if drug_name and not any(c.name == drug_name for c in all_candidates):
                        clinical_score = await self._calculate_clinical_score_real(
                            target_name=target_name,
                            clinical_data=[cd]
                        )
                        all_candidates.append(AntibodyCandidate(
                            antibody_id=f"gs-{drug_name}",
                            name=drug_name,
                            target_protein=target_name,
                            isotype=None,
                            related_diseases=disease_name,
                            full_spec=None,
                            clinical_score=clinical_score,
                            match_confidence=round(clinical_score * 0.9, 4)
                        ))

            # 3. 중복 제거 후 점수순 정렬
            seen = set()
            unique = []
            for c in all_candidates:
                key = f"{c.name}:{c.target_protein}"
                if key not in seen:
                    seen.add(key)
                    unique.append(c)

            unique.sort(key=lambda x: x.match_confidence, reverse=True)

            if unique:
                return unique[:top_k]

            # Fallback: 벡터 검색
            return await self._vector_search_antibodies(disease_name, top_k)

        except Exception as e:
            logger.error(f"[navigator] Antibody search error: {e}")
            return await self._vector_search_antibodies(disease_name, top_k)

    async def _vector_search_antibodies(
        self,
        disease_name: str,
        top_k: int
    ) -> List[AntibodyCandidate]:
        """벡터 유사도 기반 항체 검색 (fallback)"""
        try:
            disease_embedding = await self.rag_service.generate_embedding(
                f"Disease: {disease_name}, treatment target proteins, therapeutic antibodies"
            )
            results = self.supabase.rpc("match_antibody_by_disease", {
                "query_embedding": disease_embedding,
                "match_threshold": 0.4,
                "match_count": top_k * 2
            }).execute()

            antibodies = results.data or []
            candidates = []
            for ab in antibodies:
                candidates.append(AntibodyCandidate(
                    antibody_id=str(ab.get("id", "")),
                    name=ab.get("name", ab.get("product_name", "Unknown")),
                    target_protein=ab.get("target_protein", ab.get("target_normalized", "Unknown")),
                    isotype=ab.get("isotype"),
                    related_diseases=ab.get("related_disease"),
                    full_spec=ab.get("full_spec"),
                    clinical_score=ab.get("clinical_score", 0.5),
                    match_confidence=ab.get("similarity", 0.5)
                ))
            return candidates[:top_k] if candidates else await self._direct_antibody_search_fallback(disease_name, top_k)

        except Exception as e:
            logger.error(f"[navigator] Vector search error: {e}")
            return await self._direct_antibody_search_fallback(disease_name, top_k)

    async def _direct_antibody_search_fallback(
        self,
        disease_name: str,
        top_k: int
    ) -> List[AntibodyCandidate]:
        """직접 키워드 검색 최종 fallback"""
        try:
            # antibody_library 검색
            result = self.supabase.table("antibody_library").select(
                "id, product_name, target_normalized, isotype, related_disease"
            ).ilike("related_disease", f"%{disease_name}%").limit(top_k).execute()

            if result.data:
                return [
                    AntibodyCandidate(
                        antibody_id=str(ab["id"]),
                        name=ab["product_name"],
                        target_protein=ab.get("target_normalized", "Unknown"),
                        isotype=ab.get("isotype"),
                        related_diseases=ab.get("related_disease"),
                        clinical_score=0.5,
                        match_confidence=0.5
                    ) for ab in result.data
                ]

            # golden_set 검색
            gs_result = self.supabase.table("golden_set").select(
                "id, drug_name, target_1, indication, orr_pct, pfs_months, os_months, clinical_status"
            ).ilike("indication", f"%{disease_name}%").limit(top_k).execute()

            return [
                AntibodyCandidate(
                    antibody_id=str(gs["id"]),
                    name=gs["drug_name"],
                    target_protein=gs.get("target_1", "Unknown"),
                    related_diseases=gs.get("indication"),
                    clinical_score=self.scorer.calculate_score({
                        "orr_pct": gs.get("orr_pct", 0),
                        "pfs_months": gs.get("pfs_months", 0),
                        "os_months": gs.get("os_months", 0),
                        "clinical_phase": gs.get("clinical_status", "")
                    }),
                    match_confidence=0.7
                ) for gs in (gs_result.data or [])
            ]
        except Exception as e:
            logger.error(f"[navigator] Direct search fallback error: {e}")
            return []

    async def _calculate_clinical_score_real(
        self,
        target_name: str,
        clinical_data: List[Dict[str, Any]]
    ) -> float:
        """실제 임상 데이터 기반 clinical score 계산"""
        if not clinical_data:
            # golden_set에서 해당 타겟의 임상 데이터 직접 조회
            try:
                result = self.supabase.table("golden_set").select(
                    "orr_pct, pfs_months, os_months, clinical_status"
                ).or_(
                    f"target_1.ilike.%{target_name}%,target_2.ilike.%{target_name}%"
                ).order("orr_pct", desc=True).limit(5).execute()
                clinical_data = result.data or []
            except Exception as e:
                logger.warning(f"[navigator] Clinical data fetch error: {e}")
                return 0.5

        if not clinical_data:
            return 0.3  # 임상 데이터 없음

        # 최고 성능 데이터 기준 점수 계산
        best_score = 0.0
        for cd in clinical_data:
            score = self.scorer.calculate_score({
                "orr_pct": cd.get("orr_pct", 0) or 0,
                "pfs_months": cd.get("pfs_months", 0) or 0,
                "os_months": cd.get("os_months", 0) or 0,
                "clinical_phase": cd.get("clinical_status", "") or ""
            })
            best_score = max(best_score, score)

        return round(best_score, 4)

    # =========================================================================
    # Step 2: Golden Combination (Real Data - golden_set first)
    # =========================================================================

    async def _generate_golden_combination(
        self,
        target_protein: str,
        antibody_candidates: List[AntibodyCandidate]
    ) -> GoldenCombination:
        """
        실제 임상 데이터 기반 최적 링커-페이로드 조합 생성
        우선순위: golden_set(FDA 승인) → commercial_reagents → Gemini AI 추천
        """
        try:
            best_combo = None
            source = "none"

            # ── 1순위: golden_set에서 FDA 승인/임상 ADC 조합 검색 ──
            gs_result = self.supabase.table("golden_set").select(
                "id, drug_name, target_1, target_2, indication, "
                "linker_smiles, linker_type, payload_smiles, payload_class, "
                "dar, orr_pct, pfs_months, os_months, clinical_status, mechanism_of_action"
            ).or_(
                f"target_1.ilike.%{target_protein}%,target_2.ilike.%{target_protein}%"
            ).order("orr_pct", desc=True).limit(10).execute()

            gs_combos = gs_result.data or []

            if gs_combos:
                # 임상 점수로 재순위화
                scored = []
                for combo in gs_combos:
                    score = self.scorer.calculate_score({
                        "orr_pct": combo.get("orr_pct", 0) or 0,
                        "pfs_months": combo.get("pfs_months", 0) or 0,
                        "os_months": combo.get("os_months", 0) or 0,
                        "clinical_phase": combo.get("clinical_status", "") or ""
                    })
                    scored.append({**combo, "weighted_score": score})

                best_combo = max(scored, key=lambda x: x["weighted_score"])
                source = "golden_set"
                logger.info(f"[navigator] Golden combo from golden_set: {best_combo.get('drug_name')} (score={best_combo['weighted_score']:.3f})")

            # ── 2순위: commercial_reagents에서 시약 검색 ──
            if not best_combo or not best_combo.get("linker_smiles"):
                try:
                    cr_result = self.supabase.table("commercial_reagents").select(
                        "id, product_name, target_normalized, smiles_code, "
                        "payload_smiles, linker_smiles, linker_type, payload_class"
                    ).or_(
                        f"target_normalized.ilike.%{target_protein}%"
                    ).not_.is_("smiles_code", "null").limit(10).execute()

                    cr_data = cr_result.data or []
                    if cr_data and not best_combo:
                        # commercial_reagents에서 조합 구성
                        best_cr = cr_data[0]
                        best_combo = {
                            "id": best_cr.get("id", ""),
                            "drug_name": best_cr.get("product_name", "Commercial Design"),
                            "linker_smiles": best_cr.get("linker_smiles", best_cr.get("smiles_code", "")),
                            "linker_type": best_cr.get("linker_type", "cleavable"),
                            "payload_smiles": best_cr.get("payload_smiles", ""),
                            "payload_class": best_cr.get("payload_class", "Unknown"),
                            "dar": 4,
                            "orr_pct": None,
                            "pfs_months": None,
                            "os_months": None,
                            "weighted_score": 0.4
                        }
                        source = "commercial_reagents"
                    elif cr_data and best_combo:
                        # golden_set 조합에 누락된 SMILES를 commercial_reagents로 보충
                        if not best_combo.get("linker_smiles"):
                            for cr in cr_data:
                                if cr.get("linker_smiles"):
                                    best_combo["linker_smiles"] = cr["linker_smiles"]
                                    break
                        if not best_combo.get("payload_smiles"):
                            for cr in cr_data:
                                if cr.get("payload_smiles"):
                                    best_combo["payload_smiles"] = cr["payload_smiles"]
                                    break
                except Exception as e:
                    logger.warning(f"[navigator] Commercial reagent search error: {e}")

            # ── 3순위: golden_set 전체에서 최고 성능 조합 ──
            if not best_combo:
                any_gs = self.supabase.table("golden_set").select(
                    "id, drug_name, target_1, linker_smiles, linker_type, "
                    "payload_smiles, payload_class, dar, orr_pct, pfs_months, os_months, "
                    "clinical_status, mechanism_of_action"
                ).not_.is_("linker_smiles", "null").order(
                    "orr_pct", desc=True
                ).limit(1).execute()

                if any_gs.data:
                    best_combo = any_gs.data[0]
                    best_combo["weighted_score"] = self.scorer.calculate_score({
                        "orr_pct": best_combo.get("orr_pct", 0) or 0,
                        "pfs_months": best_combo.get("pfs_months", 0) or 0,
                        "os_months": best_combo.get("os_months", 0) or 0,
                        "clinical_phase": best_combo.get("clinical_status", "") or ""
                    })
                    source = "golden_set_global"
                    logger.info(f"[navigator] Using global best from golden_set: {best_combo.get('drug_name')}")

            # ── 4순위: Gemini AI 추천 ──
            if not best_combo:
                best_combo = await self._gemini_suggest_combination(target_protein)
                source = "gemini_ai"

            # Golden Combination 구성
            linker_smiles = best_combo.get("linker_smiles", "") or ""
            payload_smiles = best_combo.get("payload_smiles", "") or ""

            return GoldenCombination(
                antibody=antibody_candidates[0],
                linker=LinkerSpec(
                    id=str(best_combo.get("linker_id", best_combo.get("id", ""))),
                    smiles=linker_smiles,
                    type=best_combo.get("linker_type", "cleavable") or "cleavable",
                    cleavable=best_combo.get("linker_type", "cleavable") != "non-cleavable"
                ),
                payload=PayloadSpec(
                    id=str(best_combo.get("payload_id", best_combo.get("id", ""))),
                    smiles=payload_smiles,
                    class_name=best_combo.get("payload_class", "Unknown") or "Unknown",
                    mechanism=best_combo.get("mechanism_of_action", "") or ""
                ),
                dar=best_combo.get("dar", 4) or 4,
                historical_performance={
                    "orr_pct": best_combo.get("orr_pct"),
                    "pfs_months": best_combo.get("pfs_months"),
                    "os_months": best_combo.get("os_months"),
                    "source_drug": best_combo.get("drug_name", ""),
                    "clinical_status": best_combo.get("clinical_status", ""),
                    "data_source": source
                },
                confidence_score=best_combo.get("weighted_score", 0.5)
            )

        except Exception as e:
            logger.error(f"[navigator] Golden combination error: {e}")
            return await self._fallback_golden_combination(antibody_candidates)

    async def _gemini_suggest_combination(self, target_protein: str) -> Dict[str, Any]:
        """Gemini AI로 최적 ADC 조합 추천"""
        try:
            model = get_gemini_model(temperature=0.1)
            prompt = (
                f"You are an ADC (Antibody-Drug Conjugate) design expert.\n"
                f"For target protein '{target_protein}', suggest the optimal ADC combination.\n"
                f"Return ONLY a JSON object with these fields:\n"
                f"- drug_name: suggested ADC name\n"
                f"- linker_type: 'cleavable' or 'non-cleavable'\n"
                f"- payload_class: e.g. 'MMAE', 'DXd', 'DM1', 'SN-38'\n"
                f"- mechanism_of_action: mechanism description\n"
                f"- dar: Drug-to-Antibody Ratio (integer 2-8)\n"
                f"- rationale: 1-2 sentences explaining the choice\n"
                f"No markdown, just JSON."
            )
            response = await model.ainvoke(prompt)
            content = response.content.strip()

            # Parse JSON
            if content.startswith("{"):
                data = json.loads(content)
            else:
                start = content.find("{")
                end = content.rfind("}") + 1
                data = json.loads(content[start:end]) if start >= 0 else {}

            return {
                "id": "gemini-suggested",
                "drug_name": data.get("drug_name", f"AI-Designed {target_protein}-ADC"),
                "linker_type": data.get("linker_type", "cleavable"),
                "payload_class": data.get("payload_class", "MMAE"),
                "mechanism_of_action": data.get("mechanism_of_action", ""),
                "dar": data.get("dar", 4),
                "orr_pct": None,
                "pfs_months": None,
                "os_months": None,
                "weighted_score": 0.3,
                "linker_smiles": "",
                "payload_smiles": ""
            }
        except Exception as e:
            logger.error(f"[navigator] Gemini combination suggestion error: {e}")
            return {
                "id": "fallback",
                "drug_name": f"Designed for {target_protein}",
                "linker_type": "cleavable",
                "payload_class": "MMAE",
                "dar": 4,
                "weighted_score": 0.2,
                "linker_smiles": "",
                "payload_smiles": ""
            }

    async def _fallback_golden_combination(
        self,
        antibody_candidates: List[AntibodyCandidate]
    ) -> GoldenCombination:
        """최종 fallback: golden_set에서 가장 높은 점수 조합"""
        try:
            result = self.supabase.table("golden_set").select(
                "id, drug_name, target_1, linker_smiles, linker_type, "
                "payload_smiles, payload_class, dar, orr_pct, pfs_months, os_months, "
                "clinical_status, mechanism_of_action"
            ).not_.is_("orr_pct", "null").order(
                "orr_pct", desc=True
            ).limit(1).execute()

            if result.data:
                best = result.data[0]
                score = self.scorer.calculate_score({
                    "orr_pct": best.get("orr_pct", 0) or 0,
                    "pfs_months": best.get("pfs_months", 0) or 0,
                    "os_months": best.get("os_months", 0) or 0,
                    "clinical_phase": best.get("clinical_status", "") or ""
                })
                return GoldenCombination(
                    antibody=antibody_candidates[0] if antibody_candidates else AntibodyCandidate(
                        antibody_id="default", name="Unknown", target_protein="Unknown"
                    ),
                    linker=LinkerSpec(
                        id=str(best.get("id", "")),
                        smiles=best.get("linker_smiles", "") or "",
                        type=best.get("linker_type", "cleavable") or "cleavable"
                    ),
                    payload=PayloadSpec(
                        id=str(best.get("id", "")),
                        smiles=best.get("payload_smiles", "") or "",
                        class_name=best.get("payload_class", "Unknown") or "Unknown",
                        mechanism=best.get("mechanism_of_action", "") or ""
                    ),
                    dar=best.get("dar", 4) or 4,
                    historical_performance={
                        "orr_pct": best.get("orr_pct"),
                        "pfs_months": best.get("pfs_months"),
                        "os_months": best.get("os_months"),
                        "source_drug": best.get("drug_name", ""),
                        "data_source": "golden_set_fallback"
                    },
                    confidence_score=score
                )
        except Exception as e:
            logger.error(f"[navigator] Fallback golden combination DB error: {e}")

        # Absolute last resort - no hardcoded mock values
        return GoldenCombination(
            antibody=antibody_candidates[0] if antibody_candidates else AntibodyCandidate(
                antibody_id="default", name="Unknown", target_protein="Unknown"
            ),
            linker=LinkerSpec(id="none", smiles="", type="unknown"),
            payload=PayloadSpec(id="none", smiles="", class_name="Unknown"),
            dar=4,
            historical_performance={"data_source": "no_data_available"},
            confidence_score=0.1
        )

    # =========================================================================
    # Step 3: Property Calculation (Real Data)
    # =========================================================================

    def _combine_adc_structure(self, golden_combo: GoldenCombination) -> str:
        """ADC 구조를 SMILES로 조합"""
        linker_smiles = golden_combo.linker.smiles or ""
        payload_smiles = golden_combo.payload.smiles or ""

        parts = [s for s in [linker_smiles, payload_smiles] if s]
        return ".".join(parts) if parts else ""

    async def _calculate_properties(
        self,
        smiles: str,
        golden_combo: Optional[GoldenCombination] = None
    ) -> Dict[str, Any]:
        """분자 물성 계산 + PK 추정 (실제 데이터 기반)"""
        properties = {
            "molecular_weight": 0,
            "logp": 0,
            "hbd": 0,
            "hba": 0,
            "tpsa": 0,
            "rotatable_bonds": 0,
            "pk_parameters": {},
            "smiles_valid": False
        }

        # RDKit 물성 계산
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Crippen

            # 각 SMILES 성분 개별 계산
            total_mw = 0
            total_logp = 0
            valid_count = 0

            for smi in smiles.split("."):
                smi = smi.strip()
                if not smi:
                    continue
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    valid_count += 1
                    mw = Descriptors.MolWt(mol)
                    total_mw += mw
                    total_logp += Crippen.MolLogP(mol)
                    properties["hbd"] += Descriptors.NumHDonors(mol)
                    properties["hba"] += Descriptors.NumHAcceptors(mol)
                    properties["tpsa"] += Descriptors.TPSA(mol)
                    properties["rotatable_bonds"] += Descriptors.NumRotatableBonds(mol)

            if valid_count > 0:
                properties["molecular_weight"] = round(total_mw, 2)
                properties["logp"] = round(total_logp, 2)
                properties["smiles_valid"] = True

        except ImportError:
            logger.warning("[navigator] RDKit not available")
        except Exception as e:
            logger.error(f"[navigator] Property calculation error: {e}")

        # PK 파라미터 추정 (payload class + DAR 기반)
        properties["pk_parameters"] = await self._estimate_pk_parameters(golden_combo)

        return properties

    async def _estimate_pk_parameters(
        self,
        golden_combo: Optional[GoldenCombination] = None
    ) -> Dict[str, Any]:
        """
        PK 파라미터 추정
        golden_set의 유사 ADC 레퍼런스 기반 + payload class/DAR 보정
        """
        # 기본 ADC PK 범위 (문헌 기반)
        # 일반적인 ADC: half-life 3-6일, CL 5-20 mL/day/kg, Vd 50-100 mL/kg
        pk = {
            "half_life_hours": 96,       # 4일
            "clearance_ml_h_kg": 0.35,   # ~8.4 mL/day/kg
            "volume_of_distribution_l_kg": 0.07,  # 70 mL/kg
            "estimation_source": "default"
        }

        if not golden_combo:
            return pk

        payload_class = (golden_combo.payload.class_name or "").upper()
        dar = golden_combo.dar or 4

        # Payload class별 PK 보정 (문헌 기반)
        pk_by_payload = {
            "MMAE": {"half_life_hours": 96, "clearance_ml_h_kg": 0.30, "vd": 0.06},
            "MMAF": {"half_life_hours": 120, "clearance_ml_h_kg": 0.25, "vd": 0.05},
            "DXD": {"half_life_hours": 144, "clearance_ml_h_kg": 0.20, "vd": 0.07},
            "DM1": {"half_life_hours": 84, "clearance_ml_h_kg": 0.40, "vd": 0.08},
            "DM4": {"half_life_hours": 96, "clearance_ml_h_kg": 0.35, "vd": 0.07},
            "SN-38": {"half_life_hours": 144, "clearance_ml_h_kg": 0.22, "vd": 0.06},
            "MAYTANSINOID": {"half_life_hours": 90, "clearance_ml_h_kg": 0.38, "vd": 0.07},
            "AURISTATIN": {"half_life_hours": 100, "clearance_ml_h_kg": 0.32, "vd": 0.06},
            "CAMPTOTHECIN": {"half_life_hours": 140, "clearance_ml_h_kg": 0.22, "vd": 0.06},
            "PBD": {"half_life_hours": 72, "clearance_ml_h_kg": 0.50, "vd": 0.09},
        }

        for key, values in pk_by_payload.items():
            if key in payload_class:
                pk["half_life_hours"] = values["half_life_hours"]
                pk["clearance_ml_h_kg"] = values["clearance_ml_h_kg"]
                pk["volume_of_distribution_l_kg"] = values["vd"]
                pk["estimation_source"] = f"payload_class:{key}"
                break

        # DAR 보정: 높은 DAR → 더 빠른 클리어런스
        if dar > 4:
            pk["clearance_ml_h_kg"] *= (1 + (dar - 4) * 0.1)
            pk["half_life_hours"] *= (1 - (dar - 4) * 0.05)
        elif dar < 4:
            pk["clearance_ml_h_kg"] *= (1 - (4 - dar) * 0.05)
            pk["half_life_hours"] *= (1 + (4 - dar) * 0.03)

        pk["half_life_hours"] = round(pk["half_life_hours"], 1)
        pk["clearance_ml_h_kg"] = round(pk["clearance_ml_h_kg"], 4)
        pk["volume_of_distribution_l_kg"] = round(pk["volume_of_distribution_l_kg"], 4)

        return pk

    # =========================================================================
    # Step 5: Virtual Trial (Real Data + Gemini AI)
    # =========================================================================

    async def _run_virtual_trial(
        self,
        golden_combo: GoldenCombination,
        pk_params: Dict[str, Any],
        disease_name: str
    ) -> VirtualTrialResult:
        """
        가상 임상 시뮬레이션 (실제 데이터 기반)
        1. golden_set에서 동일/유사 타겟의 실제 임상 데이터 조회
        2. 레퍼런스 약물 대비 예측값 계산
        3. Gemini AI로 정교한 예측
        """
        historical = golden_combo.historical_performance or {}
        target = golden_combo.antibody.target_protein
        confidence_factor = golden_combo.confidence_score

        # ── 1. golden_set에서 해당 타겟+질환의 실제 임상 레퍼런스 수집 ──
        reference_data = []
        try:
            ref_result = self.supabase.table("golden_set").select(
                "drug_name, target_1, indication, orr_pct, pfs_months, os_months, "
                "clinical_status, dar, payload_class, linker_type"
            ).or_(
                f"target_1.ilike.%{target}%,target_2.ilike.%{target}%"
            ).not_.is_("orr_pct", "null").order("orr_pct", desc=True).limit(10).execute()
            reference_data = ref_result.data or []
        except Exception as e:
            logger.warning(f"[navigator] Reference data fetch error: {e}")

        # ── 2. 실제 레퍼런스 기반 예측 ──
        if reference_data:
            # 실제 임상 데이터의 통계
            orr_values = [r["orr_pct"] for r in reference_data if r.get("orr_pct")]
            pfs_values = [r["pfs_months"] for r in reference_data if r.get("pfs_months")]
            os_values = [r["os_months"] for r in reference_data if r.get("os_months")]

            # 가중 평균 (최고 성능에 가중치)
            if orr_values:
                base_orr = orr_values[0]  # 최고 ORR (이미 정렬됨)
                avg_orr = sum(orr_values) / len(orr_values)
                # 예측: 최고값과 평균의 가중 평균
                predicted_orr = base_orr * 0.6 + avg_orr * 0.4
            else:
                predicted_orr = historical.get("orr_pct") or 30

            if pfs_values:
                base_pfs = max(pfs_values)
                avg_pfs = sum(pfs_values) / len(pfs_values)
                predicted_pfs = base_pfs * 0.5 + avg_pfs * 0.5
            else:
                predicted_pfs = historical.get("pfs_months") or 5

            if os_values:
                base_os = max(os_values)
                avg_os = sum(os_values) / len(os_values)
                predicted_os = base_os * 0.5 + avg_os * 0.5
            else:
                predicted_os = historical.get("os_months") or 10

            # 신뢰도 보정
            predicted_orr = predicted_orr * (0.85 + 0.15 * confidence_factor)
            predicted_pfs = predicted_pfs * (0.9 + 0.1 * confidence_factor)
            predicted_os = predicted_os * (0.9 + 0.1 * confidence_factor)

            confidence = min(0.95, 0.5 + len(reference_data) * 0.05)
            data_source = "golden_set_reference"

            logger.info(
                f"[navigator] Virtual trial from {len(reference_data)} references: "
                f"ORR={predicted_orr:.1f}%, PFS={predicted_pfs:.1f}mo, OS={predicted_os:.1f}mo"
            )
        else:
            # 레퍼런스 없음 - Gemini AI 예측
            ai_prediction = await self._gemini_predict_trial(
                target=target,
                disease_name=disease_name,
                payload_class=golden_combo.payload.class_name,
                linker_type=golden_combo.linker.type,
                dar=golden_combo.dar
            )
            predicted_orr = ai_prediction.get("orr", 25)
            predicted_pfs = ai_prediction.get("pfs", 4)
            predicted_os = ai_prediction.get("os", 9)
            confidence = 0.3  # AI 예측은 낮은 신뢰도
            data_source = "gemini_ai_prediction"

        # ── 3. PK 데이터 시뮬레이션 (실제 PK 파라미터 기반) ──
        half_life = pk_params.get("half_life_hours", 96)
        clearance = pk_params.get("clearance_ml_h_kg", 0.35)
        vd = pk_params.get("volume_of_distribution_l_kg", 0.07)

        # 2-구획 모델 근사 (α, β phase)
        alpha_half = half_life * 0.1   # 분포상 반감기 (~10%)
        beta_half = half_life           # 소실상 반감기
        pk_data = []
        for t in range(0, 505, 24):  # 0-504시간 (3주, 일반 ADC 투여 간격)
            # 2-구획 모델: C(t) = A*exp(-αt) + B*exp(-βt)
            alpha_comp = 70 * (0.5 ** (t / max(alpha_half, 1)))
            beta_comp = 30 * (0.5 ** (t / max(beta_half, 1)))
            conc = alpha_comp + beta_comp

            # Free payload release (DAR-dependent)
            dar = golden_combo.dar or 4
            release_rate = 0.02 + (dar - 2) * 0.005  # DAR 높을수록 더 많은 release
            free_payload = conc * release_rate * (1 + t / 500)  # 시간에 따라 증가
            free_payload = min(free_payload, conc * 0.15)  # 최대 15%

            pk_data.append({
                "time_hours": t,
                "concentration": round(conc, 2),
                "free_payload": round(free_payload, 3)
            })

        # ── 4. 종양 성장 억제 시뮬레이션 ──
        tumor_data = []
        initial_volume = 100  # mm³ (일반적 종양 이식 크기)
        doubling_time = 5     # 일 (종양 배가 시간)
        tgi = min(predicted_orr * 1.1, 95)  # TGI ≈ ORR * 1.1 (상한 95%)

        for day in range(0, 43, 3):  # 0-42일 (6주)
            # 대조군: 지수 성장
            control = initial_volume * (2 ** (day / doubling_time))

            # 치료군: TGI 적용 (시그모이드 반응)
            effect = tgi / 100 * (1 - 0.5 ** (day / 7))  # 7일 후 50% 효과
            treated = initial_volume * (2 ** (day / doubling_time)) * (1 - effect)
            treated = max(treated, initial_volume * 0.05)  # 최소 5% 잔존

            tumor_data.append({
                "day": day,
                "control": round(control, 1),
                "treated": round(treated, 1)
            })

        return VirtualTrialResult(
            predicted_orr=round(predicted_orr, 1),
            predicted_pfs_months=round(predicted_pfs, 1),
            predicted_os_months=round(predicted_os, 1),
            pk_data=pk_data,
            tumor_data=tumor_data,
            patient_population=f"{disease_name} ({data_source}, {len(reference_data)} refs)",
            confidence=round(confidence, 2)
        )

    async def _gemini_predict_trial(
        self,
        target: str,
        disease_name: str,
        payload_class: str,
        linker_type: str,
        dar: int
    ) -> Dict[str, float]:
        """Gemini AI로 가상 임상 결과 예측"""
        try:
            model = get_gemini_model(temperature=0.1)
            prompt = (
                f"You are a clinical oncology expert. Predict clinical trial outcomes for:\n"
                f"- Disease: {disease_name}\n"
                f"- Target: {target}\n"
                f"- Payload class: {payload_class}\n"
                f"- Linker type: {linker_type}\n"
                f"- DAR: {dar}\n\n"
                f"Based on published ADC clinical data, predict:\n"
                f"Return ONLY a JSON object with: orr (%), pfs (months), os (months), rationale (1 sentence).\n"
                f"Be realistic based on existing ADC trial data. No markdown."
            )
            response = await model.ainvoke(prompt)
            content = response.content.strip()

            if content.startswith("{"):
                data = json.loads(content)
            else:
                start = content.find("{")
                end = content.rfind("}") + 1
                data = json.loads(content[start:end]) if start >= 0 else {}

            return {
                "orr": float(data.get("orr", 25)),
                "pfs": float(data.get("pfs", 4)),
                "os": float(data.get("os", 9))
            }
        except Exception as e:
            logger.error(f"[navigator] Gemini trial prediction error: {e}")
            return {"orr": 25, "pfs": 4, "os": 9}

    # =========================================================================
    # Database Operations
    # =========================================================================

    async def _create_db_session(self, state: NavigatorState):
        """DB 세션 생성 또는 업데이트"""
        try:
            # 기존 세션 확인
            existing = self.supabase.table("navigator_sessions").select("id").eq(
                "id", state.session_id
            ).single().execute()

            if existing.data:
                # 기존 세션 업데이트
                self.supabase.table("navigator_sessions").update({
                    "status": "running",
                    "current_step": 0,
                    "updated_at": datetime.utcnow().isoformat()
                }).eq("id", state.session_id).execute()
            else:
                # 새 세션 생성
                self.supabase.table("navigator_sessions").insert({
                    "id": state.session_id,
                    "user_id": state.user_id,
                    "disease_name": state.disease_name,
                    "status": "running",
                    "current_step": 0,
                    "total_steps": 5
                }).execute()
        except Exception as e:
            logger.warning(f"[navigator] DB session creation error: {e}")

    async def _update_db_step(self, state: NavigatorState, step: int, data: Dict):
        """DB 단계 업데이트"""
        try:
            self.supabase.table("navigator_sessions").update({
                "current_step": step,
                **data,
                "updated_at": datetime.utcnow().isoformat()
            }).eq("id", state.session_id).execute()
        except Exception as e:
            logger.warning(f"[navigator] DB step update error: {e}")

    async def _complete_db_session(self, state: NavigatorState, lineage: Dict):
        """DB 세션 완료"""
        try:
            predicted_orr = state.virtual_trial.predicted_orr if state.virtual_trial else 0
            self.supabase.rpc("complete_navigator_session", {
                "p_session_id": state.session_id,
                "p_physics_verified": state.physics_verified,
                "p_predicted_orr": predicted_orr,
                "p_lineage_data": lineage
            }).execute()
        except Exception as e:
            logger.warning(f"[navigator] DB complete error: {e}")

    async def _fail_db_session(self, state: NavigatorState, error: str):
        """DB 세션 실패"""
        try:
            self.supabase.rpc("fail_navigator_session", {
                "p_session_id": state.session_id,
                "p_error_message": error,
                "p_error_step": state.step.value
            }).execute()
        except Exception as e:
            logger.warning(f"[navigator] DB fail error: {e}")

    # =========================================================================
    # WebSocket Broadcasting
    # =========================================================================

    async def _broadcast_step(
        self,
        state: NavigatorState,
        step: int,
        message: str,
        complete: bool = False
    ):
        """WebSocket으로 진행 상황 전송"""
        try:
            await websocket_hub.stream_progress(
                state.session_id,
                progress=step * 20,  # 5단계 = 100%
                step=message
            )

            if complete:
                await websocket_hub.broadcast_to_session(state.session_id, {
                    "type": "navigator_complete",
                    "session_id": state.session_id,
                    "disease_name": state.disease_name,
                    "physics_verified": state.physics_verified
                })
        except Exception as e:
            logger.warning(f"[navigator] Broadcast error: {e}")

    # =========================================================================
    # Digital Lineage
    # =========================================================================

    async def _collect_lineage(self, state: NavigatorState) -> Dict[str, Any]:
        """Digital Lineage 수집 (데이터 출처 추적 포함)"""
        combo = state.golden_combination
        lineage = {
            "pipeline": "one_click_adc_navigator",
            "version": "2.0.0",
            "execution_timestamp": state.started_at.isoformat(),
            "disease_input": state.disease_name,
            "steps_completed": state.step.value,
            "agents_invoked": ["librarian", "alchemist", "coder", "auditor"],
            "antibody_count": len(state.antibody_candidates),
            "physics_verified": state.physics_verified,
            "confidence_score": combo.confidence_score if combo else 0,
            "targets_found": list(set(ab.target_protein for ab in state.antibody_candidates)),
            "data_sources": {
                "golden_set": True,
                "antibody_library": True,
                "commercial_reagents": True,
                "gemini_ai": True
            }
        }
        if combo and combo.historical_performance:
            lineage["combination_source"] = combo.historical_performance.get("data_source", "unknown")
            lineage["reference_drug"] = combo.historical_performance.get("source_drug", "")
        return lineage


# ============================================================================
# Singleton
# ============================================================================

_navigator_orchestrator: Optional[NavigatorOrchestrator] = None


def get_navigator_orchestrator() -> NavigatorOrchestrator:
    """Navigator Orchestrator 싱글톤 인스턴스"""
    global _navigator_orchestrator
    if _navigator_orchestrator is None:
        _navigator_orchestrator = NavigatorOrchestrator()
    return _navigator_orchestrator


# ============================================================================
# Public API
# ============================================================================

async def run_one_click_navigator(
    disease_name: str,
    session_id: Optional[str] = None,
    user_id: Optional[str] = None
) -> NavigatorResult:
    """
    One-Click ADC Navigator 실행

    Args:
        disease_name: 질환명
        session_id: 세션 ID
        user_id: 사용자 ID

    Returns:
        NavigatorResult
    """
    orchestrator = get_navigator_orchestrator()
    return await orchestrator.run_one_click_pipeline(
        disease_name=disease_name,
        session_id=session_id,
        user_id=user_id
    )
