"""
ADC Design Engine State Schema
Multi-Agent 설계 엔진용 공유 상태 스키마 (LangGraph TypedDict)
"""
from typing import TypedDict, Literal, List, Optional, Dict, Any
from datetime import datetime


class CalculatedMetrics(TypedDict, total=False):
    """계산된 물성 메트릭"""
    mw: Optional[float]           # Molecular Weight
    logp: Optional[float]         # LogP
    hbd: Optional[int]            # H-bond Donors
    hba: Optional[int]            # H-bond Acceptors
    psa: Optional[float]          # Polar Surface Area
    tpsa: Optional[float]         # Topological PSA
    rotatable_bonds: Optional[int]
    sa_score: Optional[float]     # Synthetic Accessibility


class ValidationFlags(TypedDict, total=False):
    """검증 플래그"""
    lipinski_pass: bool
    pains_free: bool
    structural_alerts: List[str]
    golden_set_similarity: float
    auditor_approved: bool


class CandidateStructure(TypedDict):
    """후보 구조"""
    smiles: str
    rank: int
    score: float
    metrics: Optional[CalculatedMetrics]
    validation: Optional[ValidationFlags]
    is_masked: bool  # Free tier 마스킹 여부


class ConstraintCheck(TypedDict):
    """제약조건 검사 결과"""
    passed: bool
    violations: List[Dict[str, Any]]
    warnings: List[Dict[str, Any]]
    constraints_checked: List[str]
    summary: str


class DesignSessionState(TypedDict):
    """
    ADC Design Engine 공유 상태

    Constraint 1 준수: 모든 에이전트가 실시간 동기화하는 상태 객체
    - current_smiles: 현재 설계 중인 SMILES
    - calculated_metrics: 계산된 MW, LogP 등
    - validation_flags: Lipinski, PAINS 등 검증 결과
    """
    # Session Info
    session_id: str
    session_type: Literal["denovo", "optimization", "audit", "cmc"]
    tier: Literal["free", "premium"]
    user_id: str

    # Design Parameters (Constraint 4: Guardrail용 초기 조건)
    target_antigen: str
    target_indication: str
    requested_dar: int
    linker_preference: str
    design_goal: str

    # Real-time Synced Data
    current_smiles: str
    calculated_metrics: CalculatedMetrics
    validation_flags: ValidationFlags

    # Workflow Control
    step: int
    total_steps: int
    current_agent: str
    status: Literal["pending", "running", "completed", "failed", "manual_review"]
    error: Optional[str]

    # Self-healing State
    requires_healing: bool
    healing_attempts: int
    last_code: Optional[str]
    last_error: Optional[str]

    # Results
    candidates: List[CandidateStructure]
    constraint_check: Optional[ConstraintCheck]
    final_report: Optional[Dict[str, Any]]

    # Librarian References
    scaffold_mappings: List[Dict[str, Any]]
    pmid_references: List[str]
    golden_set_references: List[str]

    # AlphaFold (Premium only)
    alphafold_job_id: Optional[str]
    alphafold_result: Optional[Dict[str, Any]]


def create_initial_state(
    session_id: str,
    user_id: str,
    session_type: str,
    tier: str,
    target_antigen: str = "",
    target_indication: str = "",
    requested_dar: int = 4,
    linker_preference: str = "any",
    design_goal: str = ""
) -> DesignSessionState:
    """초기 상태 생성"""
    return DesignSessionState(
        # Session Info
        session_id=session_id,
        session_type=session_type,
        tier=tier,
        user_id=user_id,

        # Design Parameters
        target_antigen=target_antigen,
        target_indication=target_indication,
        requested_dar=requested_dar,
        linker_preference=linker_preference,
        design_goal=design_goal,

        # Real-time Data
        current_smiles="",
        calculated_metrics={},
        validation_flags={},

        # Workflow
        step=0,
        total_steps=5,
        current_agent="orchestrator",
        status="pending",
        error=None,

        # Self-healing
        requires_healing=False,
        healing_attempts=0,
        last_code=None,
        last_error=None,

        # Results
        candidates=[],
        constraint_check=None,
        final_report=None,

        # References
        scaffold_mappings=[],
        pmid_references=[],
        golden_set_references=[],

        # AlphaFold
        alphafold_job_id=None,
        alphafold_result=None
    )
