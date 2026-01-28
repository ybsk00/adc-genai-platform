"""
WebSocket Hub - 실시간 상태 브로드캐스트
Constraint 1: Shared State 실시간 동기화
"""
from typing import Dict, Set, Optional
from fastapi import WebSocket
from datetime import datetime
import asyncio
import json
import logging

logger = logging.getLogger(__name__)


class WebSocketHub:
    """
    WebSocket 연결 관리 및 브로드캐스트

    세션별 구독자 관리로 설계 워크플로우 상태를 실시간 전송
    """

    def __init__(self):
        # session_id -> Set[WebSocket]
        self._connections: Dict[str, Set[WebSocket]] = {}
        self._lock = asyncio.Lock()

    async def connect(self, websocket: WebSocket, session_id: str) -> None:
        """
        WebSocket 연결 등록

        Args:
            websocket: WebSocket 연결
            session_id: 구독할 세션 ID
        """
        await websocket.accept()

        async with self._lock:
            if session_id not in self._connections:
                self._connections[session_id] = set()
            self._connections[session_id].add(websocket)

        logger.info(f"[ws] Client connected to session {session_id}")

    async def disconnect(self, websocket: WebSocket, session_id: str) -> None:
        """
        WebSocket 연결 해제

        Args:
            websocket: WebSocket 연결
            session_id: 세션 ID
        """
        async with self._lock:
            if session_id in self._connections:
                self._connections[session_id].discard(websocket)
                if not self._connections[session_id]:
                    del self._connections[session_id]

        logger.info(f"[ws] Client disconnected from session {session_id}")

    async def broadcast(self, session_id: str, message: dict) -> None:
        """
        세션 구독자들에게 메시지 브로드캐스트

        Args:
            session_id: 대상 세션 ID
            message: 전송할 메시지 (dict)
        """
        async with self._lock:
            connections = self._connections.get(session_id, set()).copy()

        if not connections:
            return

        # 타임스탬프 추가
        message["server_timestamp"] = datetime.utcnow().isoformat()

        payload = json.dumps(message, default=str)
        dead_connections = []

        for websocket in connections:
            try:
                await websocket.send_text(payload)
            except Exception as e:
                logger.warning(f"[ws] Failed to send message: {e}")
                dead_connections.append(websocket)

        # 죽은 연결 정리
        if dead_connections:
            async with self._lock:
                for ws in dead_connections:
                    if session_id in self._connections:
                        self._connections[session_id].discard(ws)

    async def broadcast_agent_status(
        self,
        session_id: str,
        agent: str,
        status: str,
        message: str,
        step: int = 0,
        extra_data: Optional[dict] = None
    ) -> None:
        """
        에이전트 상태 브로드캐스트

        Args:
            session_id: 세션 ID
            agent: 에이전트 이름
            status: 상태 (started, completed, error, warning)
            message: 상태 메시지
            step: 현재 스텝
            extra_data: 추가 데이터
        """
        payload = {
            "type": "agent_status",
            "agent": agent,
            "status": status,
            "message": message,
            "step": step,
            "timestamp": datetime.utcnow().isoformat()
        }

        if extra_data:
            payload["data"] = extra_data

        await self.broadcast(session_id, payload)

    async def broadcast_progress(
        self,
        session_id: str,
        current_step: int,
        total_steps: int,
        current_agent: str,
        progress_percent: float
    ) -> None:
        """
        진행률 브로드캐스트

        Args:
            session_id: 세션 ID
            current_step: 현재 스텝
            total_steps: 전체 스텝
            current_agent: 현재 에이전트
            progress_percent: 진행률 (0-100)
        """
        await self.broadcast(session_id, {
            "type": "progress",
            "current_step": current_step,
            "total_steps": total_steps,
            "current_agent": current_agent,
            "progress_percent": min(100, max(0, progress_percent)),
            "timestamp": datetime.utcnow().isoformat()
        })

    async def broadcast_candidates(
        self,
        session_id: str,
        candidates: list,
        is_partial: bool = False
    ) -> None:
        """
        후보 구조 브로드캐스트

        Args:
            session_id: 세션 ID
            candidates: 후보 리스트
            is_partial: 부분 결과 여부
        """
        await self.broadcast(session_id, {
            "type": "candidates_update",
            "candidates": candidates,
            "is_partial": is_partial,
            "count": len(candidates),
            "timestamp": datetime.utcnow().isoformat()
        })

    async def broadcast_completion(
        self,
        session_id: str,
        status: str,
        final_report: Optional[dict] = None
    ) -> None:
        """
        완료 브로드캐스트

        Args:
            session_id: 세션 ID
            status: 최종 상태 (completed, failed, manual_review)
            final_report: 최종 리포트
        """
        await self.broadcast(session_id, {
            "type": "session_complete",
            "status": status,
            "final_report": final_report,
            "timestamp": datetime.utcnow().isoformat()
        })

    # ====== New UI Enhancement Events ======

    async def broadcast_shared_state_sync(
        self,
        session_id: str,
        current_smiles: str,
        calculated_metrics: dict,
        scaffolds: Optional[list] = None
    ) -> None:
        """
        Shared State 동기화 브로드캐스트

        Args:
            session_id: 세션 ID
            current_smiles: 현재 SMILES
            calculated_metrics: 계산된 지표들
            scaffolds: 추출된 스캐폴드 정보
        """
        await self.broadcast(session_id, {
            "type": "shared_state_sync",
            "current_smiles": current_smiles,
            "calculated_metrics": calculated_metrics,
            "scaffolds": scaffolds or [],
            "timestamp": datetime.utcnow().isoformat()
        })

    async def broadcast_healer_action(
        self,
        session_id: str,
        action: str,
        error_type: Optional[str] = None,
        error_message: Optional[str] = None,
        fix_explanation: Optional[str] = None,
        attempt: int = 1,
        max_attempts: int = 3,
        original_code: Optional[str] = None,
        fixed_code: Optional[str] = None
    ) -> None:
        """
        Healer 자가 치유 액션 브로드캐스트

        Args:
            session_id: 세션 ID
            action: 액션 타입 (detecting, analyzing, healing, healed, failed)
            error_type: 에러 유형
            error_message: 에러 메시지
            fix_explanation: 수정 설명
            attempt: 현재 시도 횟수
            max_attempts: 최대 시도 횟수
            original_code: 원본 코드
            fixed_code: 수정된 코드
        """
        await self.broadcast(session_id, {
            "type": "healer_action",
            "action": action,
            "error_type": error_type,
            "error_message": error_message,
            "fix_explanation": fix_explanation,
            "attempt": attempt,
            "max_attempts": max_attempts,
            "original_code": original_code,
            "fixed_code": fixed_code,
            "timestamp": datetime.utcnow().isoformat()
        })

    async def broadcast_librarian_references(
        self,
        session_id: str,
        references: list,
        golden_set_refs: list
    ) -> None:
        """
        Librarian 근거 문헌 브로드캐스트

        Args:
            session_id: 세션 ID
            references: 문헌 참조 리스트
            golden_set_refs: Golden Set 참조 리스트
        """
        await self.broadcast(session_id, {
            "type": "librarian_references",
            "references": references,
            "golden_set_refs": golden_set_refs,
            "timestamp": datetime.utcnow().isoformat()
        })

    async def broadcast_auditor_feedback(
        self,
        session_id: str,
        is_in_redesign_loop: bool,
        check_items: list,
        redesign_request: Optional[dict] = None
    ) -> None:
        """
        Auditor 피드백 브로드캐스트

        Args:
            session_id: 세션 ID
            is_in_redesign_loop: Redesign Loop 상태
            check_items: 체크 항목 리스트
            redesign_request: 재설계 요청 정보
        """
        await self.broadcast(session_id, {
            "type": "auditor_feedback",
            "is_in_redesign_loop": is_in_redesign_loop,
            "check_items": check_items,
            "redesign_request": redesign_request,
            "timestamp": datetime.utcnow().isoformat()
        })

    async def broadcast_digital_seal(
        self,
        session_id: str,
        record_hash: str,
        previous_hash: Optional[str] = None,
        chain_hash: Optional[str] = None,
        is_verified: bool = True
    ) -> None:
        """
        Digital Seal 해시 정보 브로드캐스트

        Args:
            session_id: 세션 ID
            record_hash: 현재 레코드 해시
            previous_hash: 이전 해시 (체인 링크)
            chain_hash: 체인 해시
            is_verified: 검증 여부
        """
        await self.broadcast(session_id, {
            "type": "digital_seal",
            "session_id": session_id,
            "record_hash": record_hash,
            "previous_hash": previous_hash,
            "chain_hash": chain_hash,
            "is_verified": is_verified,
            "timestamp": datetime.utcnow().isoformat()
        })

    def get_connection_count(self, session_id: str) -> int:
        """세션 연결 수 조회"""
        return len(self._connections.get(session_id, set()))

    def get_active_sessions(self) -> list:
        """활성 세션 목록 조회"""
        return list(self._connections.keys())


# 싱글톤 인스턴스
websocket_hub = WebSocketHub()
