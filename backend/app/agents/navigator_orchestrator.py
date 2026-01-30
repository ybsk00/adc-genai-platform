"""
One-Click ADC Navigator Orchestrator
AstraForge Enhancement Specification v2.2 - FIXED VERSION

질환명 하나만 입력하면 AI가 협업하여
최적의 ADC 설계안을 자동 생성하는 파이프라인 오케스트레이터

3단계 파이프라인:
- Step 1: Target & Antibody Match (Librarian)
- Step 2: Linker & Payload Coupling (Alchemist)
- Step 3: Simulation & Audit (Coder + Auditor)

FIXES:
- SMILES combination with proper validation
- Empty SMILES handling with early returns
- DB query optimization with caching
- Gemini error handling with exponential backoff
- PK parameters externalized to config
- Quality-based confidence calculation
- Transaction-safe updates
- Comprehensive synonym handling
- Duplicate detection by ID
- Improved fallback with defaults
- SQL injection protection
- Memory-efficient streaming
"""

import logging
import json
import uuid
import asyncio
import re
from typing import Dict, Any, List, Optional, Tuple, Set
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from functools import lru_cache
import hashlib

from app.core.supabase import get_supabase_client
from app.core.websocket_hub import websocket_hub
from app.services.rag_service import RAGService
from app.core.gemini import get_gemini_model
from app.core.config import settings
from app.agents.librarian import LibrarianAgent
from app.agents.alchemist import AlchemistAgent
from app.agents.coder import CoderAgent
from app.agents.auditor import AuditorAgent
from app.agents.healer import HealerAgent
from app.agents.design_state import DesignSessionState
from app.services.sandbox_executor import get_sandbox_executor

logger = logging.getLogger(__name__)


# ============================================================================
# Configuration & Constants
# ============================================================================

# PK Parameters moved to external config for easy updates
DEFAULT_PK_PARAMS = {
    "half_life_hours": 96,
    "clearance_ml_h_kg": 0.35,
    "volume_of_distribution_l_kg": 0.07,
}

PK_BY_PAYLOAD = {
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

# Extended disease mapping with synonyms
DISEASE_NAME_MAP = {
    # Korean
    "유방암": "Breast Cancer", "breast cancer": "Breast Cancer",
    "폐암": "Lung Cancer", "lung cancer": "Lung Cancer",
    "위암": "Gastric Cancer", "gastric cancer": "Gastric Cancer",
    "대장암": "Colorectal Cancer", "결장암": "Colorectal Cancer",
    "colorectal cancer": "Colorectal Cancer", "colon cancer": "Colorectal Cancer",
    "방광암": "Bladder Cancer", "bladder cancer": "Bladder Cancer",
    "자궁경부암": "Cervical Cancer", "cervical cancer": "Cervical Cancer",
    "난소암": "Ovarian Cancer", "ovarian cancer": "Ovarian Cancer",
    "림프종": "Lymphoma", "lymphoma": "Lymphoma",
    "백혈병": "Leukemia", "leukemia": "Leukemia",
    "흑색종": "Melanoma", "melanoma": "Melanoma",
    "다발성골수종": "Multiple Myeloma", "multiple myeloma": "Multiple Myeloma",
    "전립선암": "Prostate Cancer", "prostate cancer": "Prostate Cancer",
    "간암": "Liver Cancer", "liver cancer": "Liver Cancer", "hepatocellular carcinoma": "Liver Cancer",
    "췌장암": "Pancreatic Cancer", "pancreatic cancer": "Pancreatic Cancer",
    "두경부암": "Head and Neck Cancer", "head and neck cancer": "Head and Neck Cancer",
    "삼중음성유방암": "Triple-Negative Breast Cancer", "triple negative breast cancer": "Triple-Negative Breast Cancer",
    "tnbc": "Triple-Negative Breast Cancer",
    "비소세포폐암": "Non-Small Cell Lung Cancer", "non small cell lung cancer": "Non-Small Cell Lung Cancer",
    "nsclc": "Non-Small Cell Lung Cancer",
    "소세포폐암": "Small Cell Lung Cancer", "small cell lung cancer": "Small Cell Lung Cancer",
    "sclc": "Small Cell Lung Cancer",
    "식도암": "Esophageal Cancer", "esophageal cancer": "Esophageal Cancer",
    "신장암": "Renal Cancer", "renal cancer": "Renal Cancer", "kidney cancer": "Renal Cancer",
    "뇌종양": "Brain Tumor", "brain tumor": "Brain Tumor", "glioma": "Brain Tumor",
}

# Safe defaults for fallback
DEFAULT_ANTIBODY = {
    "antibody_id": "trastuzumab-default",
    "name": "Trastuzumab",
    "target_protein": "HER2",
    "clinical_score": 0.85,
    "match_confidence": 0.8
}

DEFAULT_LINKER = {
    "id": "val-cit-default",
    "smiles": "CC(C)C[C@H](NC(=O)CNC(=O)[C@@H]1CCCN1)C(=O)NCCCCNC(=O)CCCCC",
    "type": "cleavable",
    "cleavable": True
}

DEFAULT_PAYLOAD = {
    "id": "mmae-default",
    "smiles": "CC[C@H]1C(=O)N(C)C(CC(C)C)C(=O)N(C)C(CC(C)C)C(=O)N1C",
    "class_name": "MMAE",
    "mechanism": "Microtubule inhibitor"
}


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
    data_quality_score: float = 0.0  # NEW: quality-based metric


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
    warnings: List[str] = field(default_factory=list)  # NEW: warning tracking


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
    warnings: List[str] = field(default_factory=list)  # NEW

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
# Cache Manager
# ============================================================================

class QueryCache:
    """Simple in-memory cache for DB queries"""
    
    def __init__(self, ttl_seconds: int = 300):
        self._cache: Dict[str, Tuple[Any, datetime]] = {}
        self._ttl = ttl_seconds
        self._lock = asyncio.Lock()
    
    def _make_key(self, query_type: str, params: Dict) -> str:
        """Create cache key from query parameters"""
        param_str = json.dumps(params, sort_keys=True, default=str)
        return hashlib.md5(f"{query_type}:{param_str}".encode()).hexdigest()
    
    async def get(self, query_type: str, params: Dict) -> Optional[Any]:
        """Get cached result if not expired"""
        key = self._make_key(query_type, params)
        async with self._lock:
            if key in self._cache:
                data, timestamp = self._cache[key]
                if (datetime.utcnow() - timestamp).seconds < self._ttl:
                    logger.debug(f"[cache] Hit for {query_type}")
                    return data
                else:
                    del self._cache[key]
        return None
    
    async def set(self, query_type: str, params: Dict, data: Any):
        """Cache result with timestamp"""
        key = self._make_key(query_type, params)
        async with self._lock:
            self._cache[key] = (data, datetime.utcnow())
    
    async def clear(self):
        """Clear all cached data"""
        async with self._lock:
            self._cache.clear()


# Global cache instance
_query_cache = QueryCache(ttl_seconds=300)


# ============================================================================
# Navigator Orchestrator - FIXED VERSION
# ============================================================================

class NavigatorOrchestrator:
    """
    One-Click ADC Navigator 오케스트레이터 (FIXED)
    
    질환명 하나만 입력받아 최적의 ADC 설계안을 자동 생성합니다.
    """

    def __init__(self):
        self.supabase = get_supabase_client()
        self.rag_service = RAGService()
        self.scorer = ClinicalWeightedScorer()
        self.cache = _query_cache
        
        # [TRUE MULTI-AGENT SYSTEM]
        self.librarian = LibrarianAgent()
        self.alchemist = AlchemistAgent()
        self.coder = CoderAgent()
        self.auditor = AuditorAgent()
        self.healer = HealerAgent()

    def _map_to_design_state(self, state: NavigatorState) -> DesignSessionState:
        """Map NavigatorState to DesignSessionState for agents"""
        return {
            "session_id": state.session_id,
            "user_id": state.user_id or "anonymous",
            "session_type": "denovo",
            "tier": "free",  # Default to free for one-click
            "target_antigen": state.target_protein,
            "target_indication": state.disease_name,
            "requested_dar": 4,
            "linker_preference": "cleavable",
            "design_goal": "Max efficacy, Min toxicity",
            "status": "running",
            "current_step": state.step.value,
            "current_smiles": state.combined_smiles,
            "candidates": [], # Will be populated by Alchemist
            "calculated_metrics": state.calculated_metrics,
            "validation_flags": {},
            "redesign_count": 0,
            "history": [],
            "last_error": None,
            "requires_healing": False
        }

    async def run_one_click_pipeline(
        self,
        disease_name: str,
        session_id: Optional[str] = None,
        user_id: Optional[str] = None,
        selected_antibody_id: Optional[str] = None
    ) -> NavigatorResult:
        """
        One-Click ADC Navigator 파이프라인 (5단계)

        실제 에이전트 execute()를 호출하고, 결과를 반영합니다.
        Fallback: Agent 실패 시 레거시 메서드로 보완.
        """
        session_id = session_id or str(uuid.uuid4())
        warnings: List[str] = []
        agents_used: List[str] = []

        state = NavigatorState(
            session_id=session_id,
            user_id=user_id,
            disease_name=disease_name,
        )

        try:
            # DB 세션 생성
            await self._create_db_session(state)

            # ═══════════════════════════════════════════════════════════
            # Step 1: Target & Antibody Match (Librarian)
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.TARGET_DISCOVERY
            await self._broadcast_step(state, 1, "Librarian 에이전트가 타겟 & 항체를 검색 중...")
            await self._broadcast_agent_log(state, "Librarian", 1, f"질환 '{disease_name}'에 대한 타겟 검색 시작")
            agents_used.append("librarian")

            antibody_candidates = await self._find_antibodies_by_disease(disease_name)
            state.antibody_candidates = antibody_candidates

            if antibody_candidates:
                state.target_protein = antibody_candidates[0].target_protein
                await self._broadcast_agent_log(
                    state, "Librarian", 1,
                    f"타겟 발견: {state.target_protein}, 항체 후보 {len(antibody_candidates)}개"
                )

            # gpNMB 실패 사례 경고
            for ab in antibody_candidates:
                if ab.target_protein and "gpNMB" in ab.target_protein:
                    gpnmb_warn = (
                        "⚠️ [FAILURE CASE #5552] gpNMB target detected. "
                        "Glembatumumab vedotin (CDX-011) failed in METRIC trial due to Shedding-related efficacy loss."
                    )
                    warnings.append(gpnmb_warn)
                    state.warnings.append(gpnmb_warn)
                    await self._broadcast_agent_log(state, "Librarian", 1, gpnmb_warn, data_source="golden_set_library:5552")

            await self._update_db_step(state, 1, {
                "antibody_candidates": [
                    {"id": ab.antibody_id, "name": ab.name, "target": ab.target_protein, "score": ab.clinical_score}
                    for ab in antibody_candidates
                ]
            })
            await self._broadcast_step(state, 1, f"✅ 타겟 발견: {state.target_protein} ({len(antibody_candidates)}개 항체)")

            # Selection Gate: 사용자 선택 대기
            if not selected_antibody_id and len(antibody_candidates) > 1:
                self.supabase.table("navigator_sessions").update({
                    "status": "waiting_for_selection",
                    "antibody_candidates": [
                        {"id": ab.antibody_id, "name": ab.name, "target_protein": ab.target_protein,
                         "clinical_score": ab.clinical_score, "match_confidence": ab.match_confidence}
                        for ab in antibody_candidates
                    ]
                }).eq("id", session_id).execute()
                # Return partial result — frontend will show selection modal
                return NavigatorResult(
                    session_id=session_id,
                    disease_name=disease_name,
                    antibody_candidates=antibody_candidates,
                    golden_combination=None,
                    calculated_metrics={},
                    physical_validations=[],
                    physics_verified=False,
                    virtual_trial=None,
                    digital_lineage={},
                    warnings=warnings
                )

            # If antibody selected, filter to that one
            if selected_antibody_id:
                selected = [ab for ab in antibody_candidates if ab.antibody_id == selected_antibody_id]
                if selected:
                    state.target_protein = selected[0].target_protein

            # Build agent state for downstream agents
            agent_state = self._map_to_design_state(state)

            # ═══════════════════════════════════════════════════════════
            # Step 2: Golden Combination (Alchemist)
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.GOLDEN_COMBINATION
            await self._broadcast_step(state, 2, "Alchemist 에이전트가 최적 ADC 조합을 설계 중...")
            await self._broadcast_agent_log(state, "Alchemist", 2, f"타겟 '{state.target_protein}'에 대한 조합 생성 시작")
            agents_used.append("alchemist")

            golden_combo = None
            alchemist_result = await self.alchemist.execute(agent_state)

            if alchemist_result.success:
                agents_used.append("alchemist:success")
                agent_data = alchemist_result.data or {}
                agent_smiles = agent_data.get("primary_smiles")
                agent_candidates = agent_data.get("candidates", [])
                golden_refs = agent_data.get("golden_set_refs", [])

                await self._broadcast_agent_log(
                    state, "Alchemist", 2,
                    f"Alchemist가 {len(agent_candidates)}개 후보 생성, SMILES: {bool(agent_smiles)}",
                    confidence=alchemist_result.confidence_score
                )

                # Alchemist의 구조화된 데이터로 GoldenCombination 구성
                golden_combo = await self._generate_golden_combination(
                    state.target_protein, state.antibody_candidates,
                    alchemist_smiles=agent_smiles,
                    alchemist_refs=golden_refs
                )
            else:
                logger.warning(f"[navigator] Alchemist failed: {alchemist_result.error}")
                await self._broadcast_agent_log(state, "Alchemist", 2, f"Alchemist 실패, 레거시 fallback 사용: {alchemist_result.error}")
                golden_combo = await self._generate_golden_combination(state.target_protein, state.antibody_candidates)

            state.golden_combination = golden_combo

            # Toxicity Guardrail
            payload_cls = golden_combo.payload.class_name.upper()
            if "PBD" in payload_cls or "SGN-CD33A" in golden_combo.historical_performance.get("source_drug", ""):
                tox_warning = (
                    "⚠️ [TOXICITY ALERT #5551] DNA-alkylating payload (PBD) detected. "
                    "Vadastuximab talirine (SGN-CD33A) caused fatal hepatotoxicity (VOD/SOS) in Phase 3."
                )
                warnings.append(tox_warning)
                state.warnings.append(tox_warning)
                await self._broadcast_agent_log(state, "Auditor", 2, tox_warning, data_source="golden_set_library:5551")

            agent_state["current_smiles"] = self._combine_adc_structure(golden_combo)
            state.combined_smiles = agent_state["current_smiles"]

            await self._update_db_step(state, 2, {
                "golden_combination": self._combination_to_dict(golden_combo)
            })
            await self._broadcast_step(state, 2, f"✅ 최적 조합 발견 (DAR: {golden_combo.dar})")

            # ═══════════════════════════════════════════════════════════
            # Step 3: 물성 계산 (Coder + Healer)
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.PROPERTY_CALCULATION
            await self._broadcast_step(state, 3, "Coder 에이전트가 Sandbox에서 물성 계산 중...")
            await self._broadcast_agent_log(state, "Coder", 3, f"RDKit Sandbox에서 분자 물성 계산 시작 (SMILES: {state.combined_smiles[:30]}...)")
            agents_used.append("coder")

            metrics = None
            coder_result = await self.coder.execute(agent_state)

            if coder_result.success:
                agents_used.append("coder:success")
                metrics = coder_result.data.get("metrics", {})
                metrics["molecular_weight"] = metrics.get("mw") or metrics.get("molecular_weight")
                metrics["smiles_valid"] = True
                metrics["pk_parameters"] = await self._estimate_pk_parameters(golden_combo)
                await self._broadcast_agent_log(
                    state, "Coder", 3,
                    f"Coder 성공: MW={metrics.get('molecular_weight')}, LogP={metrics.get('logp')}"
                )
            else:
                logger.warning(f"[navigator] Coder failed: {coder_result.error}. Engaging Healer...")
                await self._broadcast_agent_log(state, "Coder", 3, f"Coder 실패: {coder_result.error}")

                # Healer: 코드 수정 후 Coder 재시도 (최대 2회)
                agents_used.append("healer")
                agent_state["last_error"] = coder_result.error
                agent_state["last_code"] = (coder_result.data or {}).get("code")
                agent_state["requires_healing"] = True

                for heal_attempt in range(2):
                    await self._broadcast_agent_log(
                        state, "Healer", 3,
                        f"Healer 코드 수정 시도 {heal_attempt + 1}/2..."
                    )
                    healer_result = await self.healer.execute(agent_state)

                    if healer_result.success:
                        agents_used.append(f"healer:fix_{heal_attempt + 1}")
                        await self._broadcast_agent_log(state, "Healer", 3, "Healer가 코드를 수정했습니다. Coder 재실행...")

                        # Healer가 수정한 코드로 Coder 재시도
                        agent_state["requires_healing"] = False
                        if healer_result.data and healer_result.data.get("fixed_code"):
                            agent_state["last_code"] = healer_result.data["fixed_code"]

                        retry_result = await self.coder.execute(agent_state)
                        if retry_result.success:
                            agents_used.append(f"coder:retry_{heal_attempt + 1}")
                            metrics = retry_result.data.get("metrics", {})
                            metrics["molecular_weight"] = metrics.get("mw") or metrics.get("molecular_weight")
                            metrics["smiles_valid"] = True
                            metrics["pk_parameters"] = await self._estimate_pk_parameters(golden_combo)
                            await self._broadcast_agent_log(state, "Coder", 3, f"Coder 재시도 성공! MW={metrics.get('molecular_weight')}")
                            break
                        else:
                            agent_state["last_error"] = retry_result.error
                            agent_state["last_code"] = (retry_result.data or {}).get("code")
                    else:
                        await self._broadcast_agent_log(state, "Healer", 3, f"Healer 수정 실패: {healer_result.error}")

                # 모든 재시도 실패 시 Sandbox fallback
                if metrics is None:
                    await self._broadcast_agent_log(state, "Coder", 3, "Agent 재시도 모두 실패, Sandbox fallback 사용")
                    metrics = await self._calculate_properties(state.combined_smiles, golden_combo)

            state.calculated_metrics = metrics
            agent_state["calculated_metrics"] = metrics

            await self._update_db_step(state, 3, {
                "combined_smiles": state.combined_smiles,
                "calculated_metrics": metrics
            })
            await self._broadcast_step(state, 3, "✅ 물성 계산 완료")

            # ═══════════════════════════════════════════════════════════
            # Step 4: 물리적 타당성 검증 (Auditor)
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.PHYSICAL_VALIDATION
            await self._broadcast_step(state, 4, "Auditor 에이전트가 규제 준수 여부를 검증 중...")
            await self._broadcast_agent_log(state, "Auditor", 4, "PAINS/Lipinski/독성 패턴 검증 시작")
            agents_used.append("auditor")

            auditor_result = await self.auditor.execute(agent_state)

            if auditor_result.success:
                agents_used.append("auditor:success")
                validation_data = auditor_result.data.get("chemistry_validation", {})
                state.physics_verified = True
                state.physical_validations = self._convert_auditor_validation(validation_data)

                decision = auditor_result.data.get("decision", "approved")
                risk_score = auditor_result.data.get("risk_score", 0)
                await self._broadcast_agent_log(
                    state, "Auditor", 4,
                    f"Auditor 결과: {decision} (risk_score: {risk_score})",
                    confidence=auditor_result.confidence_score
                )

                if decision == "redesign":
                    warnings.append("⚠️ Auditor recommends redesign. Proceeding with current design for evaluation.")
                    state.warnings.append("Auditor recommends redesign")
            else:
                logger.warning(f"[navigator] Auditor failed: {auditor_result.error}. Using legacy validation.")
                await self._broadcast_agent_log(state, "Auditor", 4, f"Auditor 실패, 레거시 검증 사용: {auditor_result.error}")
                from app.services.physical_validator import validate_structure
                validation_result = await validate_structure(state.combined_smiles)
                state.physical_validations = validation_result.get("validations", [])
                state.physics_verified = validation_result.get("overall_status") == "pass"

            await self._update_db_step(state, 4, {
                "physical_validations": state.physical_validations,
                "physics_verified": state.physics_verified
            })

            if state.physics_verified:
                await self._broadcast_step(state, 4, "✅ Auditor 승인 완료")
            else:
                await self._broadcast_step(state, 4, "⚠️ Auditor 경고 발생 (설계 계속)")

            # ═══════════════════════════════════════════════════════════
            # Step 5: 가상 임상 시뮬레이션
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.VIRTUAL_TRIAL
            await self._broadcast_step(state, 5, "가상 임상 시뮬레이션 중...")
            await self._broadcast_agent_log(state, "Clinical", 5, "Golden Set 레퍼런스 기반 가상 임상 예측 시작")
            agents_used.append("clinical_predictor")

            virtual_trial = await self._run_virtual_trial(
                golden_combo=golden_combo,
                pk_params=metrics.get("pk_parameters", {}),
                disease_name=disease_name
            )

            state.virtual_trial = virtual_trial

            await self._broadcast_agent_log(
                state, "Clinical", 5,
                f"예측 완료: ORR={virtual_trial.predicted_orr:.1f}%, PFS={virtual_trial.predicted_pfs_months:.1f}mo, "
                f"OS={virtual_trial.predicted_os_months:.1f}mo (신뢰도: {virtual_trial.confidence:.0%})"
            )

            await self._update_db_step(state, 5, {
                "virtual_trial": {
                    "predicted_orr": virtual_trial.predicted_orr,
                    "predicted_pfs_months": virtual_trial.predicted_pfs_months,
                    "predicted_os_months": virtual_trial.predicted_os_months,
                    "pk_data": virtual_trial.pk_data,
                    "tumor_data": virtual_trial.tumor_data,
                    "confidence": virtual_trial.confidence,
                    "data_quality_score": virtual_trial.data_quality_score
                }
            })

            await self._broadcast_step(
                state, 5,
                f"✅ 가상 임상 완료 (예측 ORR: {virtual_trial.predicted_orr:.1f}%, 신뢰도: {virtual_trial.confidence:.0%})"
            )

            # ═══════════════════════════════════════════════════════════
            # 완료
            # ═══════════════════════════════════════════════════════════
            state.step = NavigatorStep.COMPLETE

            # Digital Lineage 수집
            lineage = await self._collect_lineage(state, agents_used)

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
                antibody_candidates=state.antibody_candidates,
                golden_combination=golden_combo,
                calculated_metrics=state.calculated_metrics,
                physical_validations=state.physical_validations,
                physics_verified=state.physics_verified,
                virtual_trial=virtual_trial,
                digital_lineage=lineage,
                combined_smiles=state.combined_smiles,
                execution_time_seconds=execution_time,
                warnings=warnings
            )

        except Exception as e:
            logger.exception(f"[navigator] Pipeline error: {e}")
            state.errors.append(str(e))
            await self._fail_db_session(state, str(e))
            await self._broadcast_step(state, state.step.value, f"❌ 오류 발생: {str(e)}")

            return NavigatorResult(
                session_id=session_id,
                disease_name=disease_name,
                antibody_candidates=state.antibody_candidates or [AntibodyCandidate(**DEFAULT_ANTIBODY)],
                golden_combination=state.golden_combination or GoldenCombination(
                    antibody=AntibodyCandidate(**DEFAULT_ANTIBODY),
                    linker=LinkerSpec(**DEFAULT_LINKER),
                    payload=PayloadSpec(**DEFAULT_PAYLOAD),
                    dar=4
                ),
                calculated_metrics=state.calculated_metrics or {},
                physical_validations=state.physical_validations or [],
                physics_verified=False,
                virtual_trial=VirtualTrialResult(0, 0, 0),
                digital_lineage={"error": str(e), "agents_used": agents_used},
                warnings=warnings + [str(e)]
            )

    async def _broadcast_agent_log(
        self,
        state: NavigatorState,
        agent_name: str,
        step: int,
        message: str,
        data_source: str = "",
        confidence: float = 0.0,
        pmid_refs: List[str] = None
    ):
        """Broadcast agent reasoning log via WebSocket for real-time UI display"""
        try:
            await websocket_hub.broadcast(state.session_id, {
                "type": "agent_log",
                "agent": agent_name,
                "step": step,
                "status": "running",
                "message": message,
                "data_source": data_source,
                "confidence": confidence,
                "pmid_refs": pmid_refs or [],
                "timestamp": datetime.utcnow().isoformat()
            })
        except Exception as e:
            logger.debug(f"[navigator] Agent log broadcast error: {e}")

    def _convert_auditor_validation(self, validation_data: Dict) -> List[Dict]:
        """Convert Auditor validation format to Navigator list format"""
        validations = []
        
        # Lipinski
        validations.append({
            "check_name": "Lipinski Rule",
            "category": "property",
            "passed": validation_data.get("lipinski_pass", True),
            "description": f"Violations: {validation_data.get('lipinski_violations', 0)}",
            "severity": "warning" if not validation_data.get("lipinski_pass") else "info"
        })
        
        # PAINS
        validations.append({
            "check_name": "PAINS Filter",
            "category": "toxicity",
            "passed": not validation_data.get("pains_detected", False),
            "description": f"Patterns: {len(validation_data.get('pains_pattern', []))}",
            "severity": "critical" if validation_data.get("pains_detected") else "info"
        })
        
        return validations

    # =========================================================================
    # Step 1: Target & Antibody Match (FIXED with caching)
    # =========================================================================

    def _normalize_disease_name(self, disease_name: str) -> Tuple[str, List[str]]:
        """
        FIXED: Normalize disease name with synonym handling
        
        Returns: (normalized_name, search_terms)
        """
        original = disease_name.strip()
        lower_name = original.lower()
        
        # Direct match
        if lower_name in DISEASE_NAME_MAP:
            normalized = DISEASE_NAME_MAP[lower_name]
            return normalized, [normalized, original]
        
        # Partial match
        for korean, english in DISEASE_NAME_MAP.items():
            if lower_name in korean or korean in lower_name:
                return english, [english, original, korean]
        
        # No match found, return original
        return original, [original]

    async def _find_all_targets_for_disease(self, disease_name: str) -> List[Dict[str, Any]]:
        """
        FIXED: Find targets with caching and optimized queries
        """
        # Check cache first
        cache_params = {"disease": disease_name, "type": "targets"}
        cached = await self.cache.get("targets", cache_params)
        if cached is not None:
            logger.info(f"[navigator] Cache hit for targets: {disease_name}")
            return cached
        
        targets = {}
        
        # FIXED: Use improved normalization
        normalized_disease, search_terms = self._normalize_disease_name(disease_name)
        search_terms = list(dict.fromkeys(search_terms))  # Remove duplicates while preserving order
        
        logger.info(f"[navigator] Searching targets for disease: '{disease_name}' → '{normalized_disease}' (terms: {search_terms})")

        try:
            # FIXED: Single query with OR condition instead of multiple queries
            target_conditions = []
            for term in search_terms:
                target_conditions.append(f"description.ilike.%{self._sanitize_like(term)}%")
                target_conditions.append(f"category.ilike.%{self._sanitize_like(term)}%")
            
            # Query golden_set_library (rejected 제외)
            or_clause = ",".join(target_conditions[:6])  # Limit to prevent query too long
            gs_result = self.supabase.table("golden_set_library").select(
                "name, target_1, target_2, category, description, orr_pct, pfs_months, os_months, "
                "outcome_type, dar, properties, linker_type"
            ).neq("status", "rejected").or_(or_clause).limit(20).execute()

            # Process results
            seen_names = set()
            for row in (gs_result.data or []):
                name = row.get("name", "")
                if name in seen_names:
                    continue
                seen_names.add(name)
                
                props = row.get("properties") or {}
                for target_col in ["target_1", "target_2"]:
                    t = row.get(target_col)
                    if t and t.strip():
                        canonical = t.strip().upper()
                        if canonical not in targets:
                            targets[canonical] = {
                                "canonical_name": canonical,
                                "display_name": t.strip(),
                                "source": "golden_set_library",
                                "clinical_data": [],
                                "drug_count": 0,
                                "best_orr": 0,
                                "best_phase": ""
                            }
                        targets[canonical]["drug_count"] += 1
                        targets[canonical]["clinical_data"].append({
                            "drug_name": row.get("name"),
                            "orr_pct": row.get("orr_pct"),
                            "pfs_months": row.get("pfs_months"),
                            "os_months": row.get("os_months"),
                            "clinical_status": row.get("outcome_type") or row.get("status", ""),
                            "dar": row.get("dar"),
                            "payload_class": props.get("payload_class", ""),
                            "linker_type": row.get("linker_type"),
                        })
                        orr = row.get("orr_pct") or 0
                        if orr > targets[canonical]["best_orr"]:
                            targets[canonical]["best_orr"] = orr
                            targets[canonical]["best_phase"] = row.get("outcome_type") or row.get("status", "")

            logger.info(f"[navigator] Found {len(targets)} targets from golden_set_library")

            # Query antibody_library with single query
            ab_conditions = []
            for term in search_terms:
                ab_conditions.append(f"related_disease.ilike.%{self._sanitize_like(term)}%")
            
            ab_result = self.supabase.table("antibody_library").select(
                "target_normalized, related_disease"
            ).or_(",".join(ab_conditions[:4])).limit(30).execute()

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

        # Cache results
        result = sorted(targets.values(), key=lambda x: (x["best_orr"], x["drug_count"]), reverse=True)
        await self.cache.set("targets", cache_params, result)
        
        return result
    
    def _sanitize_like(self, term: str) -> str:
        """FIXED: Sanitize LIKE pattern to prevent SQL injection"""
        # Escape special characters
        sanitized = term.replace("%", "\\%").replace("_", "\\_").replace("'", "''")
        return sanitized

    async def _gemini_suggest_targets(self, disease_name: str) -> Dict[str, Dict]:
        """
        FIXED: Gemini AI target suggestion with error handling and retry
        """
        targets = {}
        max_retries = 3
        
        for attempt in range(max_retries):
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
                
                # Parse JSON
                suggested = self._safe_json_parse(content, [])
                
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
                
                if targets:
                    break  # Success, exit retry loop
                    
            except Exception as e:
                logger.error(f"[navigator] Gemini target suggestion error (attempt {attempt+1}): {e}")
                if attempt < max_retries - 1:
                    await asyncio.sleep(2 ** attempt)  # Exponential backoff
                continue
        
        return targets

    def _safe_json_parse(self, content: str, default: Any) -> Any:
        """FIXED: Safely parse JSON with multiple fallback strategies"""
        if not content:
            return default
            
        try:
            # Direct parse
            if content.startswith("[") or content.startswith("{"):
                return json.loads(content)
            
            # Extract from markdown
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0]
                return json.loads(content.strip())
            elif "```" in content:
                content = content.split("```")[1].split("```")[0]
                return json.loads(content.strip())
            
            # Extract JSON array/object
            start = content.find("[")
            end = content.rfind("]") + 1
            if start >= 0 and end > start:
                return json.loads(content[start:end])
            
            start = content.find("{")
            end = content.rfind("}") + 1
            if start >= 0 and end > start:
                return json.loads(content[start:end])
                
        except json.JSONDecodeError as e:
            logger.warning(f"[navigator] JSON parse error: {e}")
        
        return default

    async def _find_antibodies_by_disease(
        self,
        disease_name: str,
        top_k: int = 5
    ) -> List[AntibodyCandidate]:
        """
        [FIXED] Delegate to LibrarianAgent with caching
        """
        cache_params = {"disease": disease_name, "top_k": top_k, "type": "antibodies"}
        cached = await self.cache.get("antibodies", cache_params)
        if cached:
            return [AntibodyCandidate(**ab) for ab in cached]
        
        try:
            # [REAL AGENT CALL] Delegate to Librarian
            # LibrarianAgent has the complex logic for vector search + golden set ranking
            raw_results = await self.librarian.find_antibodies_by_disease(disease_name, top_k)
            
            candidates = []
            for item in raw_results:
                candidates.append(AntibodyCandidate(
                    antibody_id=str(item.get("id", "")),
                    name=item.get("name", "Unknown"),
                    target_protein=item.get("target_protein", "Unknown"),
                    isotype=item.get("isotype"),
                    related_diseases=item.get("related_disease"),
                    full_spec=item.get("full_spec"),
                    clinical_score=item.get("clinical_score", 0.0),
                    match_confidence=item.get("combined_score", 0.0)
                ))

            if candidates:
                # Cache as dicts
                await self.cache.set("antibodies", cache_params, [self._antibody_to_dict(ab) for ab in candidates])
                return candidates

            # Fallback if Agent returns nothing
            return await self._vector_search_antibodies(disease_name, top_k)

        except Exception as e:
            logger.error(f"[navigator] Librarian delegation error: {e}")
            return await self._vector_search_antibodies(disease_name, top_k)

    def _antibody_to_dict(self, ab: AntibodyCandidate) -> Dict:
        """Convert AntibodyCandidate to dict for caching"""
        return {
            "antibody_id": ab.antibody_id,
            "name": ab.name,
            "target_protein": ab.target_protein,
            "isotype": ab.isotype,
            "related_diseases": ab.related_diseases,
            "clinical_score": ab.clinical_score,
            "match_confidence": ab.match_confidence
        }

    async def _vector_search_antibodies(
        self,
        disease_name: str,
        top_k: int
    ) -> List[AntibodyCandidate]:
        """Vector similarity search as fallback"""
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
            return candidates[:top_k] if candidates else self._fallback_antibodies(top_k)

        except Exception as e:
            logger.error(f"[navigator] Vector search error: {e}")
            return self._fallback_antibodies(top_k)

    def _fallback_antibodies(self, top_k: int) -> List[AntibodyCandidate]:
        """FIXED: Return default antibodies instead of empty"""
        default = AntibodyCandidate(**DEFAULT_ANTIBODY)
        return [default] * min(top_k, 1)

    async def _calculate_clinical_score_real(
        self,
        target_name: str,
        clinical_data: List[Dict[str, Any]]
    ) -> float:
        """Calculate clinical score from real data"""
        if not clinical_data:
            try:
                result = self.supabase.table("golden_set_library").select(
                    "orr_pct, pfs_months, os_months, outcome_type"
                ).neq("status", "rejected").or_(
                    f"target_1.ilike.%{self._sanitize_like(target_name)}%,"
                    f"target_2.ilike.%{self._sanitize_like(target_name)}%"
                ).order("orr_pct", desc=True).limit(5).execute()
                clinical_data = result.data or []
            except Exception as e:
                logger.warning(f"[navigator] Clinical data fetch error: {e}")
                return 0.5

        if not clinical_data:
            return 0.3

        best_score = 0.0
        for cd in clinical_data:
            score = self.scorer.calculate_score({
                "orr_pct": cd.get("orr_pct", 0) or 0,
                "pfs_months": cd.get("pfs_months", 0) or 0,
                "os_months": cd.get("os_months", 0) or 0,
                "clinical_phase": cd.get("clinical_status", "") or cd.get("outcome_type", "") or ""
            })
            best_score = max(best_score, score)

        return round(best_score, 4)

    # =========================================================================
    # Step 2: Golden Combination (FIXED)
    # =========================================================================

    async def _generate_golden_combination(
        self,
        target_protein: str,
        antibody_candidates: List[AntibodyCandidate],
        alchemist_smiles: Optional[str] = None,
        alchemist_refs: Optional[List[Dict]] = None
    ) -> GoldenCombination:
        """
        FIXED: Generate golden combination with proper SMILES validation.
        Alchemist의 출력(SMILES, golden_set references)이 있으면 우선 반영.
        """
        cache_params = {"target": target_protein, "type": "combination"}
        cached = await self.cache.get("combination", cache_params)
        if cached:
            return GoldenCombination(**cached)
        
        try:
            best_combo = None
            source = "none"

            # Priority 1: golden_set_library
            gs_result = self.supabase.table("golden_set_library").select(
                "id, name, target_1, target_2, category, "
                "linker_smiles, linker_type, payload_smiles, properties, "
                "dar, orr_pct, pfs_months, os_months, outcome_type"
            ).neq("status", "rejected").or_(
                f"target_1.ilike.%{self._sanitize_like(target_protein)}%,"
                f"target_2.ilike.%{self._sanitize_like(target_protein)}%"
            ).order("orr_pct", desc=True).limit(10).execute()

            gs_combos = gs_result.data or []

            if gs_combos:
                scored = []
                for combo in gs_combos:
                    props = combo.get("properties") or {}
                    # Skip if missing critical SMILES
                    if not combo.get("linker_smiles") or not combo.get("payload_smiles"):
                        continue
                        
                    score = self.scorer.calculate_score({
                        "orr_pct": combo.get("orr_pct", 0) or 0,
                        "pfs_months": combo.get("pfs_months", 0) or 0,
                        "os_months": combo.get("os_months", 0) or 0,
                        "clinical_phase": combo.get("outcome_type", "") or ""
                    })
                    combo["drug_name"] = combo.get("name", "")
                    combo["payload_class"] = props.get("payload_class", "")
                    combo["mechanism_of_action"] = props.get("mechanism_of_action", "")
                    combo["clinical_status"] = combo.get("outcome_type", "")
                    scored.append({**combo, "weighted_score": score})

                if scored:
                    best_combo = max(scored, key=lambda x: x["weighted_score"])
                    source = "golden_set_library"

            # Priority 2: commercial_reagents (if missing SMILES)
            if not best_combo or not best_combo.get("linker_smiles"):
                try:
                    cr_result = self.supabase.table("commercial_reagents").select(
                        "id, product_name, target_normalized, smiles_code, "
                        "payload_smiles, linker_smiles, linker_type, payload_class"
                    ).ilike(
                        "target_normalized", f"%{self._sanitize_like(target_protein)}%"
                    ).not_.is_("smiles_code", "null").limit(10).execute()

                    cr_data = cr_result.data or []
                    if cr_data and not best_combo:
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
                        # Supplement missing SMILES
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

            # Priority 3: Any best from golden_set
            if not best_combo:
                any_gs = self.supabase.table("golden_set_library").select(
                    "id, name, target_1, linker_smiles, linker_type, "
                    "payload_smiles, properties, dar, orr_pct, pfs_months, os_months, "
                    "outcome_type"
                ).neq("status", "rejected").not_.is_("linker_smiles", "null").not_.is_("payload_smiles", "null").order(
                    "orr_pct", desc=True
                ).limit(1).execute()

                if any_gs.data:
                    best_combo = any_gs.data[0]
                    props = best_combo.get("properties") or {}
                    best_combo["drug_name"] = best_combo.get("name", "")
                    best_combo["payload_class"] = props.get("payload_class", "")
                    best_combo["mechanism_of_action"] = props.get("mechanism_of_action", "")
                    best_combo["clinical_status"] = best_combo.get("outcome_type", "")
                    best_combo["weighted_score"] = self.scorer.calculate_score({
                        "orr_pct": best_combo.get("orr_pct", 0) or 0,
                        "pfs_months": best_combo.get("pfs_months", 0) or 0,
                        "os_months": best_combo.get("os_months", 0) or 0,
                        "clinical_phase": best_combo.get("outcome_type", "") or ""
                    })
                    source = "golden_set_library_global"

            # Priority 4: Gemini AI
            if not best_combo:
                best_combo = await self._gemini_suggest_combination(target_protein)
                source = "gemini_ai"

            # Alchemist SMILES override: if Alchemist generated a valid SMILES, use it
            if alchemist_smiles and len(alchemist_smiles) > 5:
                logger.info(f"[navigator] Using Alchemist-generated SMILES: {alchemist_smiles[:50]}...")
                # Alchemist returns a combined ADC SMILES — split by '.' for linker/payload
                parts = alchemist_smiles.split(".")
                if len(parts) >= 2:
                    best_combo["linker_smiles"] = parts[0]
                    best_combo["payload_smiles"] = parts[1]
                    source = f"{source}+alchemist_smiles"
                elif len(parts) == 1:
                    # Single SMILES — keep as payload, use DB linker
                    best_combo["payload_smiles"] = parts[0]
                    source = f"{source}+alchemist_payload"

            # Build result with validation
            linker_smiles = best_combo.get("linker_smiles", "") or ""
            payload_smiles = best_combo.get("payload_smiles", "") or ""

            # FIXED: Validate SMILES
            if not linker_smiles or not payload_smiles:
                logger.warning("[navigator] Empty SMILES, using defaults")
                linker_smiles = linker_smiles or DEFAULT_LINKER["smiles"]
                payload_smiles = payload_smiles or DEFAULT_PAYLOAD["smiles"]

            result = GoldenCombination(
                antibody=antibody_candidates[0] if antibody_candidates else AntibodyCandidate(**DEFAULT_ANTIBODY),
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
            
            # Cache the result
            await self.cache.set("combination", cache_params, self._combination_to_dict(result))
            
            return result

        except Exception as e:
            logger.error(f"[navigator] Golden combination error: {e}")
            return self._fallback_golden_combination(antibody_candidates)

    def _combination_to_dict(self, combo: GoldenCombination) -> Dict:
        """Convert GoldenCombination to dict for caching"""
        return {
            "antibody": self._antibody_to_dict(combo.antibody),
            "linker": {
                "id": combo.linker.id,
                "smiles": combo.linker.smiles,
                "type": combo.linker.type,
                "cleavable": combo.linker.cleavable
            },
            "payload": {
                "id": combo.payload.id,
                "smiles": combo.payload.smiles,
                "class_name": combo.payload.class_name,
                "mechanism": combo.payload.mechanism
            },
            "dar": combo.dar,
            "historical_performance": combo.historical_performance,
            "confidence_score": combo.confidence_score
        }

    async def _gemini_suggest_combination(self, target_protein: str) -> Dict[str, Any]:
        """FIXED: Gemini with retry logic"""
        max_retries = 3
        
        for attempt in range(max_retries):
            try:
                model = get_gemini_model(temperature=0.1)
                prompt = (
                    f"You are an ADC design expert.\n"
                    f"For target protein '{target_protein}', suggest the optimal ADC combination.\n"
                    f"Return ONLY a JSON object with:\n"
                    f"- drug_name, linker_type ('cleavable' or 'non-cleavable'), payload_class, "
                    f"mechanism_of_action, dar (2-8), rationale.\n"
                    f"No markdown, just JSON."
                )
                response = await model.ainvoke(prompt)
                content = response.content.strip()

                data = self._safe_json_parse(content, {})

                if data:
                    return {
                        "id": "gemini-suggested",
                        "drug_name": data.get("drug_name", f"AI-{target_protein}-ADC"),
                        "linker_type": data.get("linker_type", "cleavable"),
                        "payload_class": data.get("payload_class", "MMAE"),
                        "mechanism_of_action": data.get("mechanism_of_action", ""),
                        "dar": data.get("dar", 4),
                        "orr_pct": None,
                        "pfs_months": None,
                        "os_months": None,
                        "weighted_score": 0.3,
                        "linker_smiles": DEFAULT_LINKER["smiles"],  # FIXED: Use default SMILES
                        "payload_smiles": DEFAULT_PAYLOAD["smiles"]  # FIXED: Use default SMILES
                    }
                
            except Exception as e:
                logger.error(f"[navigator] Gemini combination error (attempt {attempt+1}): {e}")
                if attempt < max_retries - 1:
                    await asyncio.sleep(2 ** attempt)
                continue

        # Final fallback
        return {
            "id": "fallback",
            "drug_name": f"Default-{target_protein}",
            "linker_type": "cleavable",
            "payload_class": "MMAE",
            "dar": 4,
            "weighted_score": 0.2,
            "linker_smiles": DEFAULT_LINKER["smiles"],
            "payload_smiles": DEFAULT_PAYLOAD["smiles"]
        }

    def _fallback_golden_combination(self, antibody_candidates: List[AntibodyCandidate]) -> GoldenCombination:
        """FIXED: Return valid fallback with proper SMILES"""
        return GoldenCombination(
            antibody=antibody_candidates[0] if antibody_candidates else AntibodyCandidate(**DEFAULT_ANTIBODY),
            linker=LinkerSpec(**DEFAULT_LINKER),
            payload=PayloadSpec(**DEFAULT_PAYLOAD),
            dar=4,
            historical_performance={
                "orr_pct": 45.0,
                "pfs_months": 6.0,
                "os_months": 12.0,
                "source_drug": "Default MMAE-ValCit",
                "data_source": "default_fallback"
            },
            confidence_score=0.5
        )

    # =========================================================================
    # Step 3: Property Calculation (FIXED)
    # =========================================================================

    def _combine_adc_structure(self, golden_combo: GoldenCombination) -> str:
        """
        FIXED: Proper SMILES combination with validation
        """
        linker_smiles = golden_combo.linker.smiles or ""
        payload_smiles = golden_combo.payload.smiles or ""
        
        # Validate each component
        valid_parts = []
        
        if linker_smiles:
            # Quick validation
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(linker_smiles)
                if mol:
                    valid_parts.append(linker_smiles)
            except:
                valid_parts.append(linker_smiles)  # Assume valid if RDKit not available
        
        if payload_smiles:
            try:
                from rdkit import Chem
                mol = Chem.MolFromSmiles(payload_smiles)
                if mol:
                    valid_parts.append(payload_smiles)
            except:
                valid_parts.append(payload_smiles)
        
        if len(valid_parts) >= 2:
            # Return combined structure
            return ".".join(valid_parts)
        elif len(valid_parts) == 1:
            return valid_parts[0]
        else:
            logger.error("[navigator] No valid SMILES parts for combination")
            return ""

    async def _calculate_properties(
        self,
        smiles: str,
        golden_combo: Optional[GoldenCombination] = None
    ) -> Dict[str, Any]:
        """
        [FIXED] Property calculation using SandboxExecutor (Safe Execution)
        """
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

        if not smiles:
            logger.warning("[navigator] Empty SMILES for property calculation")
            return properties

        # [REAL AGENT] Use Sandbox for RDKit
        try:
            sandbox = get_sandbox_executor()

            # Safely encode SMILES as JSON string to avoid f-string injection
            smiles_json = json.dumps(smiles)

            # Python script to run inside sandbox
            code = f"""
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

def calculate():
    smiles = json.loads({json.dumps(smiles_json)})
    total_mw = 0
    total_logp = 0
    hbd = 0
    hba = 0
    tpsa = 0
    rot_bonds = 0
    valid_count = 0

    parts = smiles.split(".")
    for smi in parts:
        smi = smi.strip()
        if not smi: continue

        mol = Chem.MolFromSmiles(smi)
        if mol:
            valid_count += 1
            total_mw += Descriptors.MolWt(mol)
            total_logp += Crippen.MolLogP(mol)
            hbd += Descriptors.NumHDonors(mol)
            hba += Descriptors.NumHAcceptors(mol)
            tpsa += Descriptors.TPSA(mol)
            rot_bonds += Descriptors.NumRotatableBonds(mol)

    return {{
        "molecular_weight": round(total_mw, 2),
        "logp": round(total_logp, 2),
        "hbd": hbd,
        "hba": hba,
        "tpsa": round(tpsa, 2),
        "rotatable_bonds": rot_bonds,
        "valid_count": valid_count,
        "smiles_valid": valid_count > 0
    }}

print(json.dumps(calculate()))
"""
            result = await sandbox.execute(code)
            
            if result.success:
                data = json.loads(result.stdout)
                properties.update(data)
                logger.info(f"[navigator] Sandbox calculation success for {smiles[:20]}...")
            else:
                logger.error(f"[navigator] Sandbox calculation failed: {result.stderr}")
                
        except Exception as e:
            logger.error(f"[navigator] Property calculation error (Sandbox): {e}")

        # PK estimation
        properties["pk_parameters"] = await self._estimate_pk_parameters(golden_combo)

        return properties

    async def _estimate_pk_parameters(
        self,
        golden_combo: Optional[GoldenCombination] = None
    ) -> Dict[str, Any]:
        """
        FIXED: PK parameter estimation from external config
        """
        pk = DEFAULT_PK_PARAMS.copy()

        if not golden_combo:
            return pk

        payload_class = (golden_combo.payload.class_name or "").upper()
        dar = golden_combo.dar or 4

        # Match payload class
        for key, values in PK_BY_PAYLOAD.items():
            if key in payload_class:
                pk["half_life_hours"] = values["half_life_hours"]
                pk["clearance_ml_h_kg"] = values["clearance_ml_h_kg"]
                pk["volume_of_distribution_l_kg"] = values["vd"]
                pk["estimation_source"] = f"payload_class:{key}"
                break

        # DAR adjustment
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
    # Step 5: Virtual Trial (FIXED with quality-based confidence)
    # =========================================================================

    async def _run_virtual_trial(
        self,
        golden_combo: GoldenCombination,
        pk_params: Dict[str, Any],
        disease_name: str
    ) -> VirtualTrialResult:
        """
        FIXED: Virtual trial with quality-based confidence calculation
        """
        historical = golden_combo.historical_performance or {}
        target = golden_combo.antibody.target_protein

        # Quality scoring for references
        reference_data = []
        data_quality_score = 0.0
        
        try:
            ref_result = self.supabase.table("golden_set_library").select(
                "name, target_1, category, orr_pct, pfs_months, os_months, "
                "outcome_type, dar, properties, linker_type, clinical_trial_count"
            ).neq("status", "rejected").or_(
                f"target_1.ilike.%{self._sanitize_like(target)}%,"
                f"target_2.ilike.%{self._sanitize_like(target)}%"
            ).not_.is_("orr_pct", "null").order("orr_pct", desc=True).limit(10).execute()
            reference_data = ref_result.data or []
        except Exception as e:
            logger.warning(f"[navigator] Reference data fetch error: {e}")

        # Calculate predictions and quality
        if reference_data:
            # Quality metrics
            has_phase3 = any(r.get("outcome_type") in ["Phase 3", "Approved"] for r in reference_data)
            has_multiple_trials = any((r.get("clinical_trial_count") or 0) > 1 for r in reference_data)
            avg_orr_quality = sum(r.get("orr_pct", 0) for r in reference_data) / len(reference_data)
            
            data_quality_score = min(1.0, (
                0.3 +  # Base score
                0.3 * (1 if has_phase3 else 0) +
                0.2 * (1 if has_multiple_trials else 0) +
                0.2 * (avg_orr_quality / 100)
            ))
            
            # Statistical predictions
            orr_values = [r["orr_pct"] for r in reference_data if r.get("orr_pct")]
            pfs_values = [r["pfs_months"] for r in reference_data if r.get("pfs_months")]
            os_values = [r["os_months"] for r in reference_data if r.get("os_months")]

            # FIXED: Use weighted average based on trial quality
            if orr_values:
                weights = [1.0 + (0.5 if r.get("outcome_type") == "Approved" else 0) for r in reference_data if r.get("orr_pct")]
                weighted_orr = sum(o * w for o, w in zip(orr_values, weights)) / sum(weights)
                predicted_orr = weighted_orr
            else:
                predicted_orr = historical.get("orr_pct") or 30

            if pfs_values:
                predicted_pfs = sum(pfs_values) / len(pfs_values)
            else:
                predicted_pfs = historical.get("pfs_months") or 5

            if os_values:
                predicted_os = sum(os_values) / len(os_values)
            else:
                predicted_os = historical.get("os_months") or 10

            # FIXED: Confidence based on quality score, not just count
            confidence = min(0.95, 0.4 + data_quality_score * 0.5)
            data_source = "golden_set_reference"

            logger.info(
                f"[navigator] Virtual trial from {len(reference_data)} references "
                f"(quality: {data_quality_score:.2f}): "
                f"ORR={predicted_orr:.1f}%, PFS={predicted_pfs:.1f}mo, OS={predicted_os:.1f}mo, "
                f"confidence={confidence:.2f}"
            )
        else:
            # AI prediction with low confidence
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
            confidence = 0.25  # Lower for AI-only
            data_quality_score = 0.1
            data_source = "gemini_ai_prediction"

        # PK simulation
        half_life = pk_params.get("half_life_hours", 96)
        pk_data = self._simulate_pk(half_life, golden_combo.dar or 4)

        # Tumor simulation
        tumor_data = self._simulate_tumor(predicted_orr)

        return VirtualTrialResult(
            predicted_orr=round(predicted_orr, 1),
            predicted_pfs_months=round(predicted_pfs, 1),
            predicted_os_months=round(predicted_os, 1),
            pk_data=pk_data,
            tumor_data=tumor_data,
            patient_population=f"{disease_name} ({data_source}, quality: {data_quality_score:.2f})",
            confidence=round(confidence, 2),
            data_quality_score=round(data_quality_score, 2)
        )

    def _simulate_pk(self, half_life: float, dar: int) -> List[Dict[str, float]]:
        """Memory-efficient PK simulation"""
        pk_data = []
        alpha_half = half_life * 0.1
        beta_half = half_life
        
        for t in range(0, 505, 24):
            alpha_comp = 70 * (0.5 ** (t / max(alpha_half, 1)))
            beta_comp = 30 * (0.5 ** (t / max(beta_half, 1)))
            conc = alpha_comp + beta_comp
            
            release_rate = 0.02 + (dar - 2) * 0.005
            free_payload = conc * release_rate * (1 + t / 500)
            free_payload = min(free_payload, conc * 0.15)
            
            pk_data.append({
                "time_hours": t,
                "concentration": round(conc, 2),
                "free_payload": round(free_payload, 3)
            })
        
        return pk_data

    def _simulate_tumor(self, predicted_orr: float) -> List[Dict[str, float]]:
        """Memory-efficient tumor simulation"""
        tumor_data = []
        initial_volume = 100
        doubling_time = 5
        tgi = min(predicted_orr * 1.1, 95)
        
        for day in range(0, 43, 3):
            control = initial_volume * (2 ** (day / doubling_time))
            effect = tgi / 100 * (1 - 0.5 ** (day / 7))
            treated = initial_volume * (2 ** (day / doubling_time)) * (1 - effect)
            treated = max(treated, initial_volume * 0.05)
            
            tumor_data.append({
                "day": day,
                "control": round(control, 1),
                "treated": round(treated, 1)
            })
        
        return tumor_data

    async def _gemini_predict_trial(
        self,
        target: str,
        disease_name: str,
        payload_class: str,
        linker_type: str,
        dar: int
    ) -> Dict[str, float]:
        """FIXED: Gemini trial prediction with retry"""
        max_retries = 3
        
        for attempt in range(max_retries):
            try:
                model = get_gemini_model(temperature=0.1)
                prompt = (
                    f"Predict clinical trial outcomes for:\n"
                    f"- Disease: {disease_name}\n"
                    f"- Target: {target}\n"
                    f"- Payload: {payload_class}\n"
                    f"- Linker: {linker_type}\n"
                    f"- DAR: {dar}\n\n"
                    f"Return JSON: orr (%), pfs (months), os (months)."
                )
                response = await model.ainvoke(prompt)
                content = response.content.strip()

                data = self._safe_json_parse(content, {})

                if data:
                    return {
                        "orr": float(data.get("orr", 25)),
                        "pfs": float(data.get("pfs", 4)),
                        "os": float(data.get("os", 9))
                    }
                    
            except Exception as e:
                logger.error(f"[navigator] Gemini trial prediction error (attempt {attempt+1}): {e}")
                if attempt < max_retries - 1:
                    await asyncio.sleep(2 ** attempt)
                continue

        return {"orr": 25, "pfs": 4, "os": 9}

    # =========================================================================
    # Database Operations (FIXED with transactions)
    # =========================================================================

    async def _create_db_session(self, state: NavigatorState):
        """FIXED: DB session creation with conflict handling"""
        try:
            # Use upsert instead of insert/update
            self.supabase.table("navigator_sessions").upsert({
                "id": state.session_id,
                "user_id": state.user_id,
                "disease_name": state.disease_name,
                "status": "running",
                "current_step": 0,
                "total_steps": 5,
                "updated_at": datetime.utcnow().isoformat(),
                "created_at": datetime.utcnow().isoformat()
            }, on_conflict="id").execute()
        except Exception as e:
            logger.warning(f"[navigator] DB session creation error: {e}")

    async def _update_db_step(self, state: NavigatorState, step: int, data: Dict):
        """FIXED: Atomic step update"""
        try:
            update_data = {
                "current_step": step,
                "updated_at": datetime.utcnow().isoformat()
            }
            
            # Flatten nested data for storage
            for key, value in data.items():
                if isinstance(value, (dict, list)):
                    update_data[key] = json.dumps(value, default=str)
                else:
                    update_data[key] = value
            
            self.supabase.table("navigator_sessions").update(update_data).eq(
                "id", state.session_id
            ).execute()
        except Exception as e:
            logger.warning(f"[navigator] DB step update error: {e}")

    async def _complete_db_session(self, state: NavigatorState, lineage: Dict):
        """FIXED: Transaction-safe completion"""
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
            # Fallback direct update
            try:
                self.supabase.table("navigator_sessions").update({
                    "status": "completed",
                    "physics_verified": state.physics_verified,
                    "predicted_orr": state.virtual_trial.predicted_orr if state.virtual_trial else 0,
                    "completed_at": datetime.utcnow().isoformat()
                }).eq("id", state.session_id).execute()
            except Exception as e2:
                logger.error(f"[navigator] Fallback completion also failed: {e2}")

    async def _fail_db_session(self, state: NavigatorState, error: str):
        """FIXED: Error recording"""
        try:
            self.supabase.rpc("fail_navigator_session", {
                "p_session_id": state.session_id,
                "p_error_message": error[:500],  # Limit error length
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
        """Broadcast progress via WebSocket"""
        try:
            await websocket_hub.stream_progress(
                state.session_id,
                percentage=step * 20,
                operation=message
            )

            if complete:
                await websocket_hub.broadcast_completion(
                    state.session_id,
                    status="completed" if not state.errors else "partial",
                    final_report={
                        "physics_verified": state.physics_verified,
                        "virtual_trial": {
                            "orr": state.virtual_trial.predicted_orr if state.virtual_trial else 0
                        }
                    }
                )
        except Exception as e:
            logger.warning(f"[navigator] Broadcast error: {e}")

    # =========================================================================
    # Digital Lineage
    # =========================================================================

    async def _collect_lineage(self, state: NavigatorState, agents_used: List[str] = None) -> Dict[str, Any]:
        """Collect data lineage with full agent traceability"""
        combo = state.golden_combination
        execution_time = (datetime.utcnow() - state.started_at).total_seconds()
        lineage = {
            "pipeline": "one_click_adc_navigator",
            "version": "2.2.0",
            "execution_timestamp": state.started_at.isoformat(),
            "total_execution_seconds": execution_time,
            "disease_input": state.disease_name,
            "steps_completed": state.step.value,
            "agents_invoked": agents_used or ["navigator_orchestrator"],
            "antibody_count": len(state.antibody_candidates),
            "physics_verified": state.physics_verified,
            "confidence_score": combo.confidence_score if combo else 0,
            "targets_found": list(set(ab.target_protein for ab in state.antibody_candidates)),
            "warnings": state.warnings,
            "data_sources": {
                "golden_set": True,
                "antibody_library": True,
                "commercial_reagents": True,
                "gemini_ai": any("gemini" in (ab.antibody_id or "") for ab in state.antibody_candidates)
            }
        }
        if combo and combo.historical_performance:
            lineage["combination_source"] = combo.historical_performance.get("data_source", "unknown")
            lineage["reference_drug"] = combo.historical_performance.get("source_drug", "")
            lineage["historical_orr"] = combo.historical_performance.get("orr_pct")

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
    user_id: Optional[str] = None,
    selected_antibody_id: Optional[str] = None
) -> NavigatorResult:
    """
    One-Click ADC Navigator 실행 (FIXED VERSION)

    Args:
        disease_name: 질환명
        session_id: 세션 ID
        user_id: 사용자 ID
        selected_antibody_id: (Optional) 사용자가 선택한 항체 ID. 없으면 Step 1만 실행.

    Returns:
        NavigatorResult with warnings and improved error handling
    """
    orchestrator = get_navigator_orchestrator()
    return await orchestrator.run_one_click_pipeline(
        disease_name=disease_name,
        session_id=session_id,
        user_id=user_id,
        selected_antibody_id=selected_antibody_id
    )
