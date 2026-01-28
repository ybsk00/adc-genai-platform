"""
Base Agent Class for ADC Design Engine
모든 설계 에이전트의 추상 베이스 클래스
"""
from abc import ABC, abstractmethod
from typing import Optional, Dict, Any
from datetime import datetime
from dataclasses import dataclass
import logging

from app.agents.design_state import DesignSessionState
from app.core.supabase import get_supabase_client

logger = logging.getLogger(__name__)


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


class BaseDesignAgent(ABC):
    """
    설계 에이전트 추상 베이스 클래스

    모든 에이전트는 이 클래스를 상속받아 구현
    - execute(): 메인 실행 로직
    - _log_start(): 실행 시작 로그
    - _log_complete(): 실행 완료 로그
    - _log_error(): 에러 로그
    """

    name: str = "base"

    def __init__(self):
        self.supabase = get_supabase_client()

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
        execution_time_ms: int = 0
    ):
        """
        실행 완료 로그

        Constraint 2: 추론 근거(reasoning)를 기록
        """
        try:
            self.supabase.table("agent_execution_logs").insert({
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
                "created_at": datetime.utcnow().isoformat()
            }).execute()
            logger.info(f"[{self.name}] Completed for session {session_id}")
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
