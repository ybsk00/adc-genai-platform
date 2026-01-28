"""
Design Engine Orchestrator - ADC 설계 워크플로우
6개 에이전트를 LangGraph로 오케스트레이션
"""
from typing import Dict, Any, Literal
from datetime import datetime
import asyncio
import logging

from langgraph.graph import StateGraph, END

from app.agents.design_state import (
    DesignSessionState,
    create_initial_state
)
from app.agents.alchemist import AlchemistAgent
from app.agents.coder import CoderAgent
from app.agents.healer import HealerAgent
from app.agents.auditor import AuditorAgent
from app.agents.librarian import LibrarianAgent
from app.core.supabase import get_supabase_client
from app.core.websocket_hub import websocket_hub

logger = logging.getLogger(__name__)


class DesignOrchestrator:
    """
    Design Engine Orchestrator

    워크플로우:
    1. Librarian: 관련 지식 수집
    2. Alchemist: SMILES 후보 생성
    3. Coder: 물성 계산 (Docker Sandbox)
    4. Healer: 에러 시 자가 치유 (최대 3회)
    5. Auditor: 검증 및 최종 승인
    6. (조건부) Alchemist 재설계

    Constraint 1: Shared State 실시간 동기화
    """

    def __init__(self):
        self.alchemist = AlchemistAgent()
        self.coder = CoderAgent()
        self.healer = HealerAgent()
        self.auditor = AuditorAgent()
        self.librarian = LibrarianAgent()
        self.supabase = get_supabase_client()

        self.graph = self._build_graph()

    def _build_graph(self) -> StateGraph:
        """LangGraph 워크플로우 구성"""
        workflow = StateGraph(DesignSessionState)

        # 노드 추가
        workflow.add_node("librarian", self._run_librarian)
        workflow.add_node("alchemist", self._run_alchemist)
        workflow.add_node("coder", self._run_coder)
        workflow.add_node("healer", self._run_healer)
        workflow.add_node("auditor", self._run_auditor)
        workflow.add_node("manual_review", self._manual_review)

        # Entry point
        workflow.set_entry_point("librarian")

        # Librarian -> Alchemist
        workflow.add_edge("librarian", "alchemist")

        # Alchemist -> Coder
        workflow.add_edge("alchemist", "coder")

        # Coder -> Healer (on error) or Auditor (on success)
        workflow.add_conditional_edges(
            "coder",
            self._check_coder_result,
            {
                "success": "auditor",
                "error": "healer",
                "max_retries": "manual_review"
            }
        )

        # Healer -> Coder (retry)
        workflow.add_edge("healer", "coder")

        # Auditor -> END or Alchemist (redesign) or Manual Review
        workflow.add_conditional_edges(
            "auditor",
            self._check_auditor_result,
            {
                "approved": END,
                "redesign": "alchemist",
                "manual_review": "manual_review"
            }
        )

        # Manual Review -> END
        workflow.add_edge("manual_review", END)

        return workflow.compile()

    def _check_coder_result(self, state: DesignSessionState) -> str:
        """Coder 결과 확인"""
        if state.get("requires_healing", False):
            if state.get("healing_attempts", 0) >= 3:
                return "max_retries"
            return "error"
        return "success"

    def _check_auditor_result(self, state: DesignSessionState) -> str:
        """Auditor 결과 확인"""
        if not state.get("final_report"):
            return "manual_review"

        decision = state["final_report"].get("decision", {})
        action = decision.get("action", "manual_review")

        if action == "approved" or action == "approved_with_warnings":
            return "approved"
        elif action == "redesign":
            return "redesign"
        else:
            return "manual_review"

    async def _run_librarian(self, state: DesignSessionState) -> DesignSessionState:
        """Librarian 실행"""
        await self._broadcast_status(state, "librarian", "started",
                                     "Searching knowledge base for relevant references...")

        result = await self.librarian.execute(state)

        if result.success:
            state["scaffold_mappings"] = result.data.get("scaffolds", [])
            state["pmid_references"] = result.data.get("pmid_references", [])
            state["golden_set_references"] = result.data.get("golden_set_references", [])

        state["current_agent"] = "librarian"
        state["step"] = 1

        await self._broadcast_status(state, "librarian", "completed",
                                     f"Found {len(state['pmid_references'])} literature references")
        await self._update_session_db(state)

        return state

    async def _run_alchemist(self, state: DesignSessionState) -> DesignSessionState:
        """Alchemist 실행"""
        await self._broadcast_status(state, "alchemist", "started",
                                     "Designing candidate structures based on Golden Set...")

        result = await self.alchemist.execute(state)

        if result.success:
            state["current_smiles"] = result.data.get("primary_smiles", "")
            state["candidates"] = result.data.get("candidates", [])

        state["current_agent"] = "alchemist"
        state["step"] = 2

        await self._broadcast_status(state, "alchemist", "completed",
                                     f"Generated {len(state['candidates'])} candidates")
        await self._update_session_db(state)

        return state

    async def _run_coder(self, state: DesignSessionState) -> DesignSessionState:
        """Coder 실행"""
        await self._broadcast_status(state, "coder", "started",
                                     "Calculating molecular properties...")

        result = await self.coder.execute(state)

        if result.success:
            state["calculated_metrics"] = result.data.get("metrics", {})
            state["requires_healing"] = False
        else:
            state["last_code"] = result.data.get("code", "")
            state["last_error"] = result.error
            state["requires_healing"] = True

        state["current_agent"] = "coder"
        state["step"] = 3

        status = "completed" if result.success else "error"
        message = "Properties calculated" if result.success else f"Error: {result.error[:50]}..."

        await self._broadcast_status(state, "coder", status, message)
        await self._update_session_db(state)

        return state

    async def _run_healer(self, state: DesignSessionState) -> DesignSessionState:
        """Healer 실행"""
        state["healing_attempts"] = state.get("healing_attempts", 0) + 1

        await self._broadcast_status(state, "healer", "started",
                                     f"Attempting to fix error (attempt {state['healing_attempts']}/3)...")

        result = await self.healer.execute(state)

        if result.success:
            state["last_code"] = result.data.get("fixed_code", "")
            state["requires_healing"] = False

            # Knowledge Loop: 성공 시 패턴 등록은 Coder 재실행 후
        else:
            state["requires_healing"] = True

        state["current_agent"] = "healer"

        await self._broadcast_status(state, "healer",
                                     "completed" if result.success else "error",
                                     result.reasoning or "Fix attempt completed")
        await self._update_session_db(state)

        return state

    async def _run_auditor(self, state: DesignSessionState) -> DesignSessionState:
        """Auditor 실행"""
        await self._broadcast_status(state, "auditor", "started",
                                     "Validating structure and checking constraints...")

        result = await self.auditor.execute(state)

        state["validation_flags"] = result.data.get("chemistry_validation", {})
        state["constraint_check"] = result.data.get("constraint_check")
        state["final_report"] = result.data.get("final_report")

        decision = result.data.get("decision", {})
        if decision.get("approved"):
            state["status"] = "completed"
        elif decision.get("action") == "redesign":
            state["status"] = "running"
        else:
            state["status"] = "manual_review"

        state["current_agent"] = "auditor"
        state["step"] = 4

        await self._broadcast_status(state, "auditor",
                                     "completed" if decision.get("approved") else "warning",
                                     decision.get("reasoning", "Audit completed"))
        await self._update_session_db(state)

        return state

    async def _manual_review(self, state: DesignSessionState) -> DesignSessionState:
        """Manual Review 상태로 전환"""
        state["status"] = "manual_review"
        state["current_agent"] = "manual_review"

        await self._broadcast_status(state, "system", "warning",
                                     "Session requires manual expert review")
        await self._update_session_db(state)

        return state

    async def _broadcast_status(
        self,
        state: DesignSessionState,
        agent: str,
        status: str,
        message: str
    ):
        """WebSocket으로 실시간 상태 전송"""
        logger.info(f"[{agent}] {status}: {message}")

        # WebSocket broadcast
        await websocket_hub.broadcast_agent_status(
            session_id=state["session_id"],
            agent=agent,
            status=status,
            message=message,
            step=state.get("step", 0)
        )

    async def _update_session_db(self, state: DesignSessionState):
        """세션 상태를 DB에 저장"""
        try:
            self.supabase.table("design_sessions").update({
                "status": state.get("status", "running"),
                "current_step": state.get("step", 0),
                "current_agent": state.get("current_agent"),
                "current_smiles": state.get("current_smiles"),
                "calculated_metrics": state.get("calculated_metrics", {}),
                "validation_flags": state.get("validation_flags", {}),
                "result_candidates": state.get("candidates", []),
                "final_report": state.get("final_report"),
                "updated_at": datetime.utcnow().isoformat()
            }).eq("id", state["session_id"]).execute()
        except Exception as e:
            logger.error(f"[orchestrator] Failed to update session: {e}")

    async def run(self, session_id: str) -> DesignSessionState:
        """
        설계 워크플로우 실행

        Args:
            session_id: 세션 ID (DB에서 조회)

        Returns:
            DesignSessionState: 최종 상태
        """
        # 세션 정보 로드
        session = self.supabase.table("design_sessions").select("*").eq(
            "id", session_id
        ).single().execute()

        if not session.data:
            raise ValueError(f"Session not found: {session_id}")

        # 초기 상태 생성
        initial_state = create_initial_state(
            session_id=session_id,
            user_id=session.data["user_id"],
            session_type=session.data["session_type"],
            tier=session.data["tier"],
            target_antigen=session.data.get("target_antigen", ""),
            target_indication=session.data.get("target_indication", ""),
            requested_dar=session.data.get("requested_dar", 4),
            linker_preference=session.data.get("linker_preference", "any"),
            design_goal=session.data.get("design_goal", "")
        )

        initial_state["status"] = "running"

        # 세션 상태 업데이트
        self.supabase.table("design_sessions").update({
            "status": "running"
        }).eq("id", session_id).execute()

        # 워크플로우 실행
        try:
            final_state = await self.graph.ainvoke(initial_state)

            # 완료 상태 업데이트
            self.supabase.table("design_sessions").update({
                "status": final_state.get("status", "completed"),
                "completed_at": datetime.utcnow().isoformat() if final_state.get("status") == "completed" else None
            }).eq("id", session_id).execute()

            return final_state

        except Exception as e:
            logger.exception(f"[orchestrator] Workflow error: {e}")

            # 실패 상태 업데이트
            self.supabase.table("design_sessions").update({
                "status": "failed"
            }).eq("id", session_id).execute()

            raise


async def run_design_workflow(session_id: str) -> DesignSessionState:
    """설계 워크플로우 실행 헬퍼 함수"""
    orchestrator = DesignOrchestrator()
    return await orchestrator.run(session_id)
