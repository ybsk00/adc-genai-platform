"""
Base Agent Class for ADC Design Engine
모든 설계 에이전트의 추상 베이스 클래스

Phase 2 Enhancement:
- Digital Lineage (데이터 계보) 지원
- Calculation Metadata 자동 수집
- 21 CFR Part 11 준수 강화
"""
from abc import ABC, abstractmethod
from typing import Optional, Dict, Any
from datetime import datetime
from dataclasses import dataclass, field
import logging
import hashlib
import platform
import sys
import uuid

from app.agents.design_state import DesignSessionState
from app.core.supabase import get_supabase_client
from app.core.config import settings

logger = logging.getLogger(__name__)


# ============================================================================
# Digital Lineage Metadata
# ============================================================================

@dataclass
class CalculationMetadata:
    """연산 메타데이터 (Digital Lineage)"""
    calculation_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    timestamp: str = field(default_factory=lambda: datetime.utcnow().isoformat())

    # 실행 환경
    execution_environment: Dict[str, Any] = field(default_factory=dict)

    # 라이브러리 버전
    library_versions: Dict[str, str] = field(default_factory=dict)

    # NVIDIA NIM API 정보 (Stage 3)
    nvidia_nim: Optional[Dict[str, Any]] = None

    # 시뮬레이션 파라미터 (MD용)
    simulation_params: Optional[Dict[str, Any]] = None

    # 입력/출력 해시 (재현성 검증용)
    input_hash: Optional[str] = None
    output_hash: Optional[str] = None

    # 에이전트 정보
    agent_name: str = ""
    agent_version: str = "v1.0.0"
    llm_model: Optional[str] = None
    llm_temperature: Optional[float] = None


def get_execution_environment() -> Dict[str, Any]:
    """현재 실행 환경 정보 수집"""
    return {
        "sandbox_type": "subprocess",  # TODO: Docker 지원 시 변경
        "container_image": "adc-sandbox:1.0.0",
        "python_version": sys.version.split()[0],
        "os_info": f"{platform.system()} {platform.release()}",
        "machine": platform.machine(),
        "processor": platform.processor() or "unknown"
    }


def get_library_versions() -> Dict[str, str]:
    """주요 라이브러리 버전 수집"""
    versions = {}

    # RDKit
    try:
        from rdkit import __version__ as rdkit_version
        versions["rdkit"] = rdkit_version
    except ImportError:
        versions["rdkit"] = "not_installed"

    # NumPy
    try:
        import numpy as np
        versions["numpy"] = np.__version__
    except ImportError:
        versions["numpy"] = "not_installed"

    # SciPy
    try:
        import scipy
        versions["scipy"] = scipy.__version__
    except ImportError:
        pass

    # LangChain
    try:
        import langchain
        versions["langchain"] = langchain.__version__
    except ImportError:
        pass

    # Google GenAI
    try:
        import google.generativeai as genai
        versions["google_genai"] = "installed"
    except ImportError:
        pass

    return versions


def calculate_hash(data: Any) -> str:
    """데이터 SHA-256 해시 계산"""
    import json
    if isinstance(data, dict):
        data_str = json.dumps(data, sort_keys=True, default=str)
    else:
        data_str = str(data)
    return hashlib.sha256(data_str.encode()).hexdigest()


@dataclass
class AgentOutput:
    """에이전트 실행 결과"""
    success: bool
    data: Dict[str, Any]
    reasoning: Optional[str] = None
    error: Optional[str] = None
    confidence_score: Optional[float] = None
    next_agent: Optional[str] = None
    referenced_knowledge_ids: Optional[list] = None
    referenced_pmids: Optional[list] = None

    # Phase 2: Digital Lineage
    calculation_metadata: Optional[CalculationMetadata] = None


class BaseDesignAgent(ABC):
    """
    설계 에이전트 추상 베이스 클래스

    모든 에이전트는 이 클래스를 상속받아 구현
    - execute(): 메인 실행 로직
    - _log_start(): 실행 시작 로그
    - _log_complete(): 실행 완료 로그
    - _log_error(): 에러 로그

    Phase 2 Enhancement:
    - Digital Lineage 메타데이터 자동 수집
    - 입력/출력 해시 계산
    - 실행 환경 기록
    """

    name: str = "base"
    version: str = "v1.0.0"

    def __init__(self):
        self.supabase = get_supabase_client()
        # Phase 2: 실행 환경 캐시
        self._execution_env = None
        self._library_versions = None

    def _get_execution_env(self) -> Dict[str, Any]:
        """실행 환경 정보 (캐시됨)"""
        if self._execution_env is None:
            self._execution_env = get_execution_environment()
        return self._execution_env

    def _get_library_versions(self) -> Dict[str, str]:
        """라이브러리 버전 정보 (캐시됨)"""
        if self._library_versions is None:
            self._library_versions = get_library_versions()
        return self._library_versions

    def _create_calculation_metadata(
        self,
        input_data: Optional[Dict] = None,
        output_data: Optional[Dict] = None,
        llm_model: Optional[str] = None,
        llm_temperature: Optional[float] = None,
        nvidia_nim_info: Optional[Dict] = None,
        simulation_params: Optional[Dict] = None
    ) -> CalculationMetadata:
        """Digital Lineage용 Calculation Metadata 생성"""
        return CalculationMetadata(
            calculation_id=str(uuid.uuid4()),
            timestamp=datetime.utcnow().isoformat(),
            execution_environment=self._get_execution_env(),
            library_versions=self._get_library_versions(),
            nvidia_nim=nvidia_nim_info,
            simulation_params=simulation_params,
            input_hash=calculate_hash(input_data) if input_data else None,
            output_hash=calculate_hash(output_data) if output_data else None,
            agent_name=self.name,
            agent_version=self.version,
            llm_model=llm_model,
            llm_temperature=llm_temperature
        )

    @abstractmethod
    async def execute(self, state: DesignSessionState) -> AgentOutput:
        """
        에이전트 메인 실행 로직

        Args:
            state: 현재 세션 상태

        Returns:
            AgentOutput: 실행 결과
        """
        pass

    async def _log_start(self, session_id: str, input_data: Optional[Dict] = None):
        """
        실행 시작 로그 (21 CFR Part 11 준수)

        Constraint 6: agent_execution_logs에 기록
        """
        try:
            self.supabase.table("agent_execution_logs").insert({
                "session_id": session_id,
                "agent_name": self.name,
                "status": "started",
                "input_data": input_data,
                "created_at": datetime.utcnow().isoformat()
            }).execute()
            logger.info(f"[{self.name}] Started for session {session_id}")
        except Exception as e:
            logger.error(f"[{self.name}] Failed to log start: {e}")

    async def _log_complete(
        self,
        session_id: str,
        output: AgentOutput,
        execution_time_ms: int = 0,
        input_data: Optional[Dict] = None
    ):
        """
        실행 완료 로그

        Constraint 2: 추론 근거(reasoning)를 기록

        Phase 2 Enhancement:
        - Digital Lineage 메타데이터 포함
        - input_hash/output_hash 자동 계산
        """
        try:
            # Phase 2: Calculation Metadata 생성 또는 사용
            if output.calculation_metadata:
                calc_meta = output.calculation_metadata
            else:
                calc_meta = self._create_calculation_metadata(
                    input_data=input_data,
                    output_data=output.data
                )

            log_data = {
                "session_id": session_id,
                "agent_name": self.name,
                "status": "completed",
                "reasoning": output.reasoning,
                "decision_summary": output.data.get("summary", ""),
                "confidence_score": output.confidence_score,
                "output_data": output.data,
                "execution_time_ms": execution_time_ms,
                "referenced_knowledge_ids": output.referenced_knowledge_ids or [],
                "referenced_pmids": output.referenced_pmids or [],
                "created_at": datetime.utcnow().isoformat(),
                # Phase 2: Digital Lineage
                "calculation_id": calc_meta.calculation_id,
                "agent_version": calc_meta.agent_version,
                "execution_env": calc_meta.execution_environment,
                "library_versions": calc_meta.library_versions,
                "nvidia_nim_info": calc_meta.nvidia_nim,
                "simulation_params": calc_meta.simulation_params,
                "input_hash": calc_meta.input_hash,
                "output_hash": calc_meta.output_hash,
                "llm_model": calc_meta.llm_model,
                "llm_temperature": calc_meta.llm_temperature
            }

            self.supabase.table("agent_execution_logs").insert(log_data).execute()
            logger.info(f"[{self.name}] Completed for session {session_id} (calc_id: {calc_meta.calculation_id[:8]}...)")
        except Exception as e:
            logger.error(f"[{self.name}] Failed to log complete: {e}")

    async def _log_error(
        self,
        session_id: str,
        error_message: str,
        input_data: Optional[Dict] = None
    ):
        """에러 로그"""
        try:
            self.supabase.table("agent_execution_logs").insert({
                "session_id": session_id,
                "agent_name": self.name,
                "status": "error",
                "error_message": error_message,
                "input_data": input_data,
                "created_at": datetime.utcnow().isoformat()
            }).execute()
            logger.error(f"[{self.name}] Error for session {session_id}: {error_message}")
        except Exception as e:
            logger.error(f"[{self.name}] Failed to log error: {e}")

    async def _log_healing(
        self,
        session_id: str,
        retry_count: int,
        fix_logs: list,
        healing_successful: bool
    ):
        """
        The Healer 전용: 자가 치유 로그

        Constraint 2: 수정 횟수와 에러 메시지를 기록하여 모델 튜닝 데이터로 활용
        """
        try:
            self.supabase.table("agent_execution_logs").insert({
                "session_id": session_id,
                "agent_name": "healer",
                "status": "healed" if healing_successful else "error",
                "retry_count": retry_count,
                "fix_logs": fix_logs,
                "healing_successful": healing_successful,
                "created_at": datetime.utcnow().isoformat()
            }).execute()
        except Exception as e:
            logger.error(f"[healer] Failed to log healing: {e}")
