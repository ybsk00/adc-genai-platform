"""
ADC State Schema - 에이전트 공유 메모리
모든 에이전트가 공유하는 상태 스키마 (Pydantic 모델)
"""
from typing import Optional, List, Dict, Any
from pydantic import BaseModel, Field
from enum import Enum
from datetime import datetime


class JobStatus(str, Enum):
    """시뮬레이션 작업 상태"""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    PARTIAL = "partial"  # 부분 실패


class AgentStatus(str, Enum):
    """개별 에이전트 상태"""
    PENDING = "pending"
    RUNNING = "running"
    DONE = "done"
    ERROR = "error"
    SKIPPED = "skipped"


class ToxicityRisk(BaseModel):
    """독성 리스크 정보"""
    risk_type: str = Field(..., description="리스크 유형 (e.g., Neutropenia, Hepatotoxicity)")
    severity: str = Field(..., description="심각도 (Low/Medium/High)")
    confidence: float = Field(..., ge=0, le=1, description="신뢰도 (0-1)")
    evidence: Optional[str] = None


class PatentInfo(BaseModel):
    """특허 정보"""
    patent_id: str
    title: str
    status: str  # Active, Expired, Pending
    expiry_date: Optional[str] = None
    risk_level: str = "Low"  # Low, Medium, High


class CompetitorInfo(BaseModel):
    """경쟁사 정보"""
    company: str
    drug_name: Optional[str] = None
    phase: str  # Phase 1, Phase 2, Phase 3, Approved
    target: str
    notes: Optional[str] = None


class StructureAnalysis(BaseModel):
    """구조 분석 결과"""
    stability_score: float = Field(..., ge=0, le=100)
    conjugation_sites: List[str] = Field(default_factory=list)
    aggregation_risk: str = "Low"
    pdb_data: Optional[str] = None  # PDB 구조 데이터
    analysis_notes: Optional[str] = None  # [NEW] 상세 분석 노트 (RMSD 등)



class CommercialInfo(BaseModel):
    """[NEW] 상용화 분석 정보 (Ambeed Data)"""
    feasibility_score: float = Field(..., description="구매 용이성 점수 (0-100)")
    total_estimated_cost: str = Field(..., description="초기 연구용 시약 예상 비용")
    payload_info: Optional[Dict[str, Any]] = None  # Ambeed Cat#, Price
    linker_info: Optional[Dict[str, Any]] = None
    conjugate_info: Optional[Dict[str, Any]] = None  # 링커-페이로드 접합체
    availability: str = "Unknown"  # 재고 상태



class AgentResult(BaseModel):
    """개별 에이전트 실행 결과"""
    agent_id: str
    status: AgentStatus = AgentStatus.PENDING
    started_at: Optional[datetime] = None
    completed_at: Optional[datetime] = None
    error_message: Optional[str] = None
    result: Optional[Dict[str, Any]] = None


class ADCInput(BaseModel):
    """ADC 분석 입력 데이터 (프론트엔드에서 전달)"""
    # Step 1: Target & Antibody
    target_name: str = Field(..., description="타겟 이름 (e.g., HER2, TROP-2)")
    antibody_type: str = Field(..., description="항체 타입 (e.g., trastuzumab, custom)")
    custom_sequence: Optional[str] = None  # FASTA 서열
    
    # Step 2: Payload & Linker
    payload_id: str = Field(..., description="약물 ID (e.g., mmae, dxd)")
    linker_id: str = Field(..., description="링커 ID (e.g., val-cit, mcc)")
    dar: int = Field(default=4, ge=1, le=8, description="Drug-to-Antibody Ratio")
    
    # Step 3: Configuration
    mode: str = Field(default="deep", description="분석 모드 (fast/deep)")
    job_name: Optional[str] = None


class ADCState(BaseModel):
    """
    LangGraph 공유 상태 스키마
    모든 에이전트가 이 상태를 읽고 쓸 수 있음
    """
    # 작업 메타데이터
    job_id: str
    user_id: str
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    status: JobStatus = JobStatus.PENDING
    
    # 입력 데이터
    input: ADCInput
    
    # 에이전트 상태 추적
    agents: Dict[str, AgentResult] = Field(default_factory=lambda: {
        "structure": AgentResult(agent_id="structure"),
        "toxicology": AgentResult(agent_id="toxicology"),
        "patent": AgentResult(agent_id="patent"),
        "competitor": AgentResult(agent_id="competitor"),
        "clinical": AgentResult(agent_id="clinical"),
        "report": AgentResult(agent_id="report"),
        "commercial": AgentResult(agent_id="commercial"),  # [NEW] 추가
    })

    
    # 분석 결과 (각 에이전트가 채움)
    structure_analysis: Optional[StructureAnalysis] = None
    toxicity_risks: List[ToxicityRisk] = Field(default_factory=list)
    patent_landscape: List[PatentInfo] = Field(default_factory=list)
    competitors: List[CompetitorInfo] = Field(default_factory=list)
    clinical_protocol: Optional[Dict[str, Any]] = None
    
    # [NEW] 상용화 분석 결과 (견적서)
    commercial_feasibility: Optional[CommercialInfo] = None

    
    # 최종 리포트
    final_grade: Optional[str] = None  # A, B+, B, C, D
    recommendation: Optional[str] = None  # Go, Conditional Go, No Go
    executive_summary: Optional[str] = None
    report_url: Optional[str] = None  # S3 PDF URL
    
    # 점수 (0-100)
    scores: Dict[str, float] = Field(default_factory=lambda: {
        "efficacy": 0,
        "toxicity": 0,
        "properties": 0,
        "patent": 0,
        "market": 0,
    })
    
    # 에러 추적
    errors: List[str] = Field(default_factory=list)
    
    def get_progress(self) -> int:
        """전체 진행률 계산 (0-100)"""
        agent_weights = {
            "structure": 15,
            "toxicology": 20,
            "patent": 15,
            "competitor": 15,
            "clinical": 15,
            "report": 20,
            "commercial": 10,  # [NEW] 가중치 추가 (총합 100 조정 필요)
        }

        
        completed = 0
        for agent_id, weight in agent_weights.items():
            if self.agents[agent_id].status == AgentStatus.DONE:
                completed += weight
            elif self.agents[agent_id].status == AgentStatus.RUNNING:
                completed += weight * 0.5  # 진행 중은 절반
        
        return int(completed)
    
    def update_agent_status(self, agent_id: str, status: AgentStatus, error: Optional[str] = None):
        """에이전트 상태 업데이트"""
        if agent_id in self.agents:
            self.agents[agent_id].status = status
            if status == AgentStatus.RUNNING:
                self.agents[agent_id].started_at = datetime.utcnow()
            elif status in [AgentStatus.DONE, AgentStatus.ERROR]:
                self.agents[agent_id].completed_at = datetime.utcnow()
            if error:
                self.agents[agent_id].error_message = error
                self.errors.append(f"[{agent_id}] {error}")
        self.updated_at = datetime.utcnow()
