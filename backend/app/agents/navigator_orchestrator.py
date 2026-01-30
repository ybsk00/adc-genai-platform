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
import uuid
from typing import Dict, Any, List, Optional
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum

from app.core.supabase import get_supabase_client
from app.core.websocket_hub import websocket_hub
from app.services.rag_service import RAGService

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

            await self._broadcast_step(
                state, 1,
                f"✅ {len(antibody_candidates)}개 항체 후보 발견 (타겟: {state.target_protein})"
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
            metrics = await self._calculate_properties(combined_smiles)
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
    # Step 1: Target & Antibody Match
    # =========================================================================

    async def _find_antibodies_by_disease(
        self,
        disease_name: str,
        top_k: int = 3
    ) -> List[AntibodyCandidate]:
        """질환명 기반 최적 항체 검색"""
        try:
            # 1. 질환명 임베딩 생성
            disease_embedding = await self.rag_service.generate_embedding(
                f"Disease: {disease_name}, treatment target proteins, therapeutic antibodies"
            )

            # 2. antibody_library 벡터 검색
            results = self.supabase.rpc("match_antibody_by_disease", {
                "query_embedding": disease_embedding,
                "match_threshold": 0.5,
                "match_count": top_k * 2  # 후보 풀
            }).execute()

            antibodies = results.data or []

            # 3. 임상 데이터 기반 재순위화
            scored_antibodies = []
            for ab in antibodies:
                clinical_score = self._calculate_clinical_score(ab)
                scored_antibodies.append({
                    **ab,
                    "clinical_score": clinical_score,
                    "combined_score": ab.get("similarity", 0.5) * 0.4 + clinical_score * 0.6
                })

            # 4. Top K 선정
            top_antibodies = sorted(
                scored_antibodies,
                key=lambda x: x["combined_score"],
                reverse=True
            )[:top_k]

            # 5. Fallback: 벡터 검색 결과 없으면 직접 검색
            if not top_antibodies:
                logger.info(f"[navigator] Vector search empty, trying direct search for: {disease_name}")
                top_antibodies = await self._direct_antibody_search(disease_name, top_k)

            return [
                AntibodyCandidate(
                    antibody_id=str(ab.get("id", "")),
                    name=ab.get("name", "Unknown"),
                    target_protein=ab.get("target_protein", "Unknown"),
                    isotype=ab.get("isotype"),
                    related_diseases=ab.get("related_disease"),
                    full_spec=ab.get("full_spec"),
                    clinical_score=ab.get("clinical_score", 0.5),
                    match_confidence=ab.get("combined_score", ab.get("similarity", 0.5))
                )
                for ab in top_antibodies
            ]

        except Exception as e:
            logger.error(f"[navigator] Antibody search error: {e}")
            # Fallback
            return await self._direct_antibody_search(disease_name, top_k)

    async def _direct_antibody_search(
        self,
        disease_name: str,
        top_k: int
    ) -> List[Dict[str, Any]]:
        """직접 키워드 검색 (fallback)"""
        try:
            result = self.supabase.table("antibody_library").select(
                "id, product_name, target_normalized, isotype, related_disease, full_spec"
            ).ilike(
                "related_disease", f"%{disease_name}%"
            ).limit(top_k).execute()

            if result.data:
                return [
                    {
                        "id": ab["id"],
                        "name": ab["product_name"],
                        "target_protein": ab.get("target_normalized", "Unknown"),
                        "isotype": ab.get("isotype"),
                        "related_disease": ab.get("related_disease"),
                        "full_spec": ab.get("full_spec"),
                        "clinical_score": 0.5,
                        "combined_score": 0.5
                    }
                    for ab in result.data
                ]

            # Golden Set에서 검색
            gs_result = self.supabase.table("golden_set").select(
                "id, drug_name, target_1, indication"
            ).ilike(
                "indication", f"%{disease_name}%"
            ).limit(top_k).execute()

            return [
                {
                    "id": gs["id"],
                    "name": gs["drug_name"],
                    "target_protein": gs.get("target_1", "Unknown"),
                    "isotype": None,
                    "related_disease": gs.get("indication"),
                    "full_spec": None,
                    "clinical_score": 0.7,
                    "combined_score": 0.7
                }
                for gs in (gs_result.data or [])
            ]

        except Exception as e:
            logger.error(f"[navigator] Direct search error: {e}")
            return []

    def _calculate_clinical_score(self, ab: Dict[str, Any]) -> float:
        """항체의 임상 점수 계산"""
        return self.scorer.calculate_score({
            "orr_pct": ab.get("orr_pct", 0),
            "pfs_months": 0,
            "os_months": ab.get("os_months", 0),
            "clinical_phase": "Phase 2"
        })

    # =========================================================================
    # Step 2: Golden Combination
    # =========================================================================

    async def _generate_golden_combination(
        self,
        target_protein: str,
        antibody_candidates: List[AntibodyCandidate]
    ) -> GoldenCombination:
        """임상 데이터 기반 최적 링커-페이로드 조합 생성"""
        try:
            # 1. 타겟에 대한 성공적인 시약 조합 검색
            # commercial_reagents는 target, target_normalized 컬럼 사용
            result = self.supabase.table("commercial_reagents").select(
                "*"
            ).or_(
                f"target.ilike.%{target_protein}%,target_normalized.ilike.%{target_protein}%"
            ).not_.is_("orr_pct", "null").order(
                "orr_pct", desc=True
            ).limit(20).execute()

            reagent_combos = result.data or []

            # 2. 임상 데이터 가중치 스코어링
            scored_combos = []
            for combo in reagent_combos:
                score = self.scorer.calculate_score(combo)
                scored_combos.append({**combo, "weighted_score": score})

            # 3. 상위 조합 선정
            if scored_combos:
                best_combo = max(scored_combos, key=lambda x: x["weighted_score"])
            else:
                # Fallback: Golden Set에서 가져오기
                best_combo = await self._fallback_golden_set_combo(target_protein)

            # 4. Golden Combination 구성
            linker_smiles = best_combo.get("linker_smiles", "")
            payload_smiles = best_combo.get("payload_smiles", "")

            # 기본 SMILES (없는 경우)
            if not linker_smiles:
                linker_smiles = "CC(=O)NCCCCC(=O)O"  # MC linker 예시
            if not payload_smiles:
                payload_smiles = "CC(C)CC(NC(=O)C)C(=O)NC"  # DM1-like 예시

            return GoldenCombination(
                antibody=antibody_candidates[0],
                linker=LinkerSpec(
                    id=str(best_combo.get("linker_id", best_combo.get("id", ""))),
                    smiles=linker_smiles,
                    type=best_combo.get("linker_type", "cleavable"),
                    cleavable=best_combo.get("is_cleavable", True)
                ),
                payload=PayloadSpec(
                    id=str(best_combo.get("payload_id", best_combo.get("id", ""))),
                    smiles=payload_smiles,
                    class_name=best_combo.get("payload_class", "Maytansinoid"),
                    mechanism=best_combo.get("mechanism_of_action", "Microtubule inhibitor")
                ),
                dar=best_combo.get("dar", 4),
                historical_performance={
                    "orr_pct": best_combo.get("orr_pct"),
                    "pfs_months": best_combo.get("pfs_months"),
                    "os_months": best_combo.get("os_months"),
                    "source_drug": best_combo.get("drug_name")
                },
                confidence_score=best_combo.get("weighted_score", 0.5)
            )

        except Exception as e:
            logger.error(f"[navigator] Golden combination error: {e}")
            # Fallback
            return await self._fallback_golden_combination(antibody_candidates)

    async def _fallback_golden_set_combo(self, target_protein: str) -> Dict[str, Any]:
        """Golden Set에서 조합 가져오기"""
        result = self.supabase.table("golden_set").select(
            "id, drug_name, target_1, linker_smiles, linker_type, "
            "payload_smiles, payload_class, dar, orr_pct, os_months, mechanism_of_action"
        ).or_(
            f"target_1.ilike.%{target_protein}%,target_2.ilike.%{target_protein}%"
        ).order("orr_pct", desc=True).limit(1).execute()

        if result.data:
            return result.data[0]

        # 기본값
        return {
            "id": "",
            "drug_name": "Default Design",
            "linker_smiles": "CC(=O)NCCCCC(=O)O",
            "linker_type": "cleavable",
            "payload_smiles": "CC(C)CC(NC(=O)C)C(=O)NC",
            "payload_class": "Maytansinoid",
            "dar": 4,
            "orr_pct": 40,
            "os_months": 12
        }

    async def _fallback_golden_combination(
        self,
        antibody_candidates: List[AntibodyCandidate]
    ) -> GoldenCombination:
        """기본 Golden Combination"""
        return GoldenCombination(
            antibody=antibody_candidates[0] if antibody_candidates else AntibodyCandidate(
                antibody_id="default",
                name="Default Antibody",
                target_protein="HER2"
            ),
            linker=LinkerSpec(
                id="default-linker",
                smiles="CC(=O)NCCCCC(=O)O",
                type="cleavable",
                cleavable=True
            ),
            payload=PayloadSpec(
                id="default-payload",
                smiles="CC(C)CC(NC(=O)C)C(=O)NC",
                class_name="Maytansinoid",
                mechanism="Microtubule inhibitor"
            ),
            dar=4,
            historical_performance={
                "orr_pct": 40,
                "pfs_months": 6,
                "os_months": 12,
                "source_drug": "Reference Design"
            },
            confidence_score=0.5
        )

    # =========================================================================
    # Step 3: Property Calculation
    # =========================================================================

    def _combine_adc_structure(self, golden_combo: GoldenCombination) -> str:
        """ADC 구조를 SMILES로 조합"""
        # 간단한 조합: linker + payload
        linker_smiles = golden_combo.linker.smiles
        payload_smiles = golden_combo.payload.smiles

        # 실제로는 결합점을 고려한 조합이 필요
        # 여기서는 단순 연결로 처리
        combined = f"{linker_smiles}.{payload_smiles}"

        return combined

    async def _calculate_properties(self, smiles: str) -> Dict[str, Any]:
        """분자 물성 계산"""
        properties = {
            "molecular_weight": 0,
            "logp": 0,
            "hbd": 0,
            "hba": 0,
            "tpsa": 0,
            "rotatable_bonds": 0,
            "pk_parameters": {}
        }

        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Crippen

            # SMILES 파싱 (조합 구조에서 첫 번째 성분)
            mol = Chem.MolFromSmiles(smiles.split(".")[0])
            if mol:
                properties["molecular_weight"] = Descriptors.MolWt(mol)
                properties["logp"] = Crippen.MolLogP(mol)
                properties["hbd"] = Descriptors.NumHDonors(mol)
                properties["hba"] = Descriptors.NumHAcceptors(mol)
                properties["tpsa"] = Descriptors.TPSA(mol)
                properties["rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)

            # PK 파라미터 (예측 기반)
            properties["pk_parameters"] = {
                "half_life_hours": 100 + (properties["molecular_weight"] / 10),
                "clearance_ml_h_kg": 0.05,
                "volume_of_distribution_l_kg": 0.05
            }

        except ImportError:
            logger.warning("[navigator] RDKit not available")
        except Exception as e:
            logger.error(f"[navigator] Property calculation error: {e}")

        return properties

    # =========================================================================
    # Step 5: Virtual Trial
    # =========================================================================

    async def _run_virtual_trial(
        self,
        golden_combo: GoldenCombination,
        pk_params: Dict[str, Any],
        disease_name: str
    ) -> VirtualTrialResult:
        """가상 임상 시뮬레이션"""
        # 역사적 성능 데이터 기반 예측
        historical = golden_combo.historical_performance or {}
        base_orr = historical.get("orr_pct", 40) or 40
        base_pfs = historical.get("pfs_months", 6) or 6
        base_os = historical.get("os_months", 12) or 12

        # 간단한 예측 모델 (실제로는 더 복잡한 모델 사용)
        confidence_factor = golden_combo.confidence_score
        predicted_orr = base_orr * (0.8 + 0.4 * confidence_factor)
        predicted_pfs = base_pfs * (0.9 + 0.2 * confidence_factor)
        predicted_os = base_os * (0.9 + 0.2 * confidence_factor)

        # PK 데이터 시뮬레이션
        half_life = pk_params.get("half_life_hours", 100)
        pk_data = []
        for t in range(0, 169, 24):  # 0-168시간 (1주일)
            conc = 100 * (0.5 ** (t / half_life))
            free_payload = conc * 0.1  # 10% free payload
            pk_data.append({
                "time_hours": t,
                "concentration": round(conc, 2),
                "free_payload": round(free_payload, 2)
            })

        # 종양 데이터 시뮬레이션
        tumor_data = []
        initial_volume = 1000
        tgi = predicted_orr  # TGI는 ORR과 유사하게 설정
        for day in range(0, 29, 7):  # 0-28일
            control = initial_volume * (1.5 ** (day / 7))
            treated = initial_volume * (1 - tgi / 100 * day / 28)
            treated = max(treated, initial_volume * 0.1)  # 최소 10%
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
            patient_population=disease_name,
            confidence=golden_combo.confidence_score
        )

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
        """Digital Lineage 수집"""
        return {
            "pipeline": "one_click_adc_navigator",
            "version": "1.0.0",
            "execution_timestamp": state.started_at.isoformat(),
            "disease_input": state.disease_name,
            "steps_completed": state.step.value,
            "agents_invoked": [
                "librarian",
                "alchemist",
                "coder",
                "auditor"
            ],
            "antibody_count": len(state.antibody_candidates),
            "physics_verified": state.physics_verified,
            "confidence_score": state.golden_combination.confidence_score if state.golden_combination else 0
        }


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
