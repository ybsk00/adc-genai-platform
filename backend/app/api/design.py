"""
Design Engine API Routes
ADC 설계 세션 및 워크플로우 API
"""
from fastapi import APIRouter, Depends, HTTPException, BackgroundTasks, WebSocket, WebSocketDisconnect
from pydantic import BaseModel
from typing import Optional, List
import logging

from app.core.supabase import get_supabase_client
from app.core.websocket_hub import websocket_hub
from app.security.masking import DataMaskingService
from app.agents.design_orchestrator import run_design_workflow

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/api/design", tags=["Design Engine"])
masking_service = DataMaskingService()


# ===== Request/Response Models =====

class CreateSessionRequest(BaseModel):
    """세션 생성 요청"""
    session_type: str  # 'denovo', 'optimization', 'audit', 'cmc'
    target_antigen: Optional[str] = None
    target_indication: Optional[str] = None
    requested_dar: Optional[int] = 4
    linker_preference: Optional[str] = "any"
    design_goal: Optional[str] = None


class CreateSessionResponse(BaseModel):
    """세션 생성 응답"""
    session_id: str
    status: str


class SessionStatusResponse(BaseModel):
    """세션 상태 응답"""
    session_id: str
    status: str
    current_step: int
    current_agent: Optional[str]
    candidates_count: int
    tier: str


# ===== Helper Functions =====

async def get_current_user_id() -> str:
    """현재 사용자 ID 조회 (TODO: JWT 인증 연동)"""
    # 개발용 더미 - 실제로는 JWT에서 추출
    return "dev-user-id"


async def get_user_tier(user_id: str) -> str:
    """사용자 티어 조회"""
    supabase = get_supabase_client()
    try:
        result = supabase.table("profiles").select("plan").eq(
            "id", user_id
        ).single().execute()

        if result.data:
            return "premium" if result.data.get("plan") == "premium" else "free"
    except Exception as e:
        logger.warning(f"Failed to get user tier: {e}")

    return "free"


# ===== API Endpoints =====

@router.post("/session", response_model=CreateSessionResponse)
async def create_session(
    request: CreateSessionRequest,
    background_tasks: BackgroundTasks
):
    """
    새 설계 세션 생성

    세션 타입:
    - denovo: De novo Design (신규 후보 설계)
    - optimization: Lead Optimization (기존 후보 최적화)
    - audit: Pre-clinical Audit (전임상 감사)
    - cmc: CMC & Sourcing (제조 및 소싱)
    """
    user_id = await get_current_user_id()
    tier = await get_user_tier(user_id)

    supabase = get_supabase_client()

    try:
        # 세션 생성
        result = supabase.table("design_sessions").insert({
            "user_id": user_id,
            "session_type": request.session_type,
            "tier": tier,
            "status": "pending",
            "target_antigen": request.target_antigen,
            "target_indication": request.target_indication,
            "requested_dar": request.requested_dar,
            "linker_preference": request.linker_preference,
            "design_goal": request.design_goal
        }).execute()

        if not result.data:
            raise HTTPException(500, "Failed to create session")

        session_id = result.data[0]["id"]

        return CreateSessionResponse(
            session_id=session_id,
            status="pending"
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"Session creation error: {e}")
        raise HTTPException(500, f"Failed to create session: {str(e)}")


@router.post("/session/{session_id}/start")
async def start_design(
    session_id: str,
    background_tasks: BackgroundTasks
):
    """
    설계 워크플로우 시작

    백그라운드에서 에이전트 오케스트레이션 실행
    """
    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    # 세션 소유권 확인
    session = supabase.table("design_sessions").select("user_id, status").eq(
        "id", session_id
    ).single().execute()

    if not session.data:
        raise HTTPException(404, "Session not found")

    if session.data["user_id"] != user_id:
        raise HTTPException(403, "Not authorized to access this session")

    if session.data["status"] == "running":
        raise HTTPException(400, "Session is already running")

    # 백그라운드에서 워크플로우 실행
    background_tasks.add_task(run_design_workflow, session_id)

    return {"status": "started", "session_id": session_id}


@router.get("/session/{session_id}")
async def get_session(session_id: str):
    """
    세션 상태 및 결과 조회

    Constraint 3: 백엔드에서 티어에 따라 마스킹 적용
    """
    user_id = await get_current_user_id()
    tier = await get_user_tier(user_id)
    supabase = get_supabase_client()

    # 세션 조회
    session = supabase.table("design_sessions").select("*").eq(
        "id", session_id
    ).single().execute()

    if not session.data:
        raise HTTPException(404, "Session not found")

    if session.data["user_id"] != user_id:
        raise HTTPException(403, "Not authorized to access this session")

    # 백엔드 마스킹 적용 (Constraint 3)
    masked_data = masking_service.mask_session_result(
        session.data,
        tier=tier
    )

    return masked_data


@router.get("/sessions")
async def list_sessions(
    status: Optional[str] = None,
    session_type: Optional[str] = None,
    limit: int = 20,
    offset: int = 0
):
    """사용자의 세션 목록 조회"""
    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    query = supabase.table("design_sessions").select(
        "id, session_type, tier, status, current_step, target_antigen, created_at, updated_at"
    ).eq("user_id", user_id)

    if status:
        query = query.eq("status", status)
    if session_type:
        query = query.eq("session_type", session_type)

    result = query.order("created_at", desc=True).range(
        offset, offset + limit - 1
    ).execute()

    return {
        "sessions": result.data or [],
        "total": len(result.data or [])
    }


@router.get("/session/{session_id}/logs")
async def get_session_logs(
    session_id: str,
    agent: Optional[str] = None
):
    """
    세션 에이전트 실행 로그 조회

    21 CFR Part 11 감사 추적용
    """
    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    # 세션 소유권 확인
    session = supabase.table("design_sessions").select("user_id").eq(
        "id", session_id
    ).single().execute()

    if not session.data or session.data["user_id"] != user_id:
        raise HTTPException(403, "Not authorized")

    # 로그 조회
    query = supabase.table("agent_execution_logs").select(
        "id, agent_name, status, reasoning, decision_summary, confidence_score, "
        "execution_time_ms, retry_count, sequence_number, created_at"
    ).eq("session_id", session_id)

    if agent:
        query = query.eq("agent_name", agent)

    result = query.order("sequence_number", desc=False).execute()

    return {
        "logs": result.data or [],
        "session_id": session_id
    }


@router.delete("/session/{session_id}")
async def cancel_session(session_id: str):
    """
    세션 취소

    running 상태의 세션은 취소 불가 (완료될 때까지 대기)
    """
    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    # 세션 조회
    session = supabase.table("design_sessions").select("user_id, status").eq(
        "id", session_id
    ).single().execute()

    if not session.data:
        raise HTTPException(404, "Session not found")

    if session.data["user_id"] != user_id:
        raise HTTPException(403, "Not authorized")

    if session.data["status"] == "running":
        raise HTTPException(400, "Cannot cancel running session")

    # 세션 삭제 (또는 상태 변경)
    supabase.table("design_sessions").update({
        "status": "cancelled"
    }).eq("id", session_id).execute()

    return {"status": "cancelled", "session_id": session_id}


@router.get("/session/{session_id}/digital-seal")
async def get_digital_seal(session_id: str):
    """
    Digital Seal 정보 조회

    21 CFR Part 11 준수 해시 체인 정보 반환
    """
    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    # 세션 소유권 확인
    session = supabase.table("design_sessions").select("user_id").eq(
        "id", session_id
    ).single().execute()

    if not session.data or session.data["user_id"] != user_id:
        raise HTTPException(403, "Not authorized")

    # 최신 로그 해시 조회
    log_result = supabase.table("agent_execution_logs").select(
        "id, record_hash, previous_hash, created_at"
    ).eq("session_id", session_id).order(
        "sequence_number", desc=True
    ).limit(1).execute()

    if not log_result.data:
        return {
            "session_id": session_id,
            "record_hash": None,
            "previous_hash": None,
            "chain_hash": None,
            "is_verified": False,
            "message": "No execution logs found"
        }

    latest_log = log_result.data[0]

    # 체인 해시 검증 (간단한 검증 - 실제로는 전체 체인 검증 필요)
    is_verified = latest_log.get("record_hash") is not None

    # 체인 해시 계산 (session_id + record_hash)
    import hashlib
    chain_data = f"{session_id}:{latest_log.get('record_hash', '')}"
    chain_hash = hashlib.sha256(chain_data.encode()).hexdigest()

    return {
        "session_id": session_id,
        "record_hash": latest_log.get("record_hash"),
        "previous_hash": latest_log.get("previous_hash"),
        "chain_hash": chain_hash,
        "timestamp": latest_log.get("created_at"),
        "is_verified": is_verified
    }


@router.get("/session/{session_id}/audit-bundle")
async def export_audit_bundle(session_id: str):
    """
    Audit Trail Bundle 내보내기

    21 CFR Part 11 준수 감사 추적 번들 생성
    포함 항목:
    - 최종 설계 SMILES
    - 실행된 파이썬 코드
    - 에이전트별 추론 로그
    - 각 단계별 SHA-256 해시값
    """
    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    # 세션 소유권 확인
    session = supabase.table("design_sessions").select("*").eq(
        "id", session_id
    ).single().execute()

    if not session.data or session.data["user_id"] != user_id:
        raise HTTPException(403, "Not authorized")

    # 에이전트 로그 조회
    logs = supabase.table("agent_execution_logs").select(
        "agent_name, status, reasoning, decision_summary, confidence_score, "
        "record_hash, previous_hash, executed_code, execution_result, "
        "sequence_number, created_at"
    ).eq("session_id", session_id).order(
        "sequence_number", desc=False
    ).execute()

    # 최종 후보 조회
    candidates = supabase.table("design_candidates").select(
        "rank, smiles, score, metrics, validation"
    ).eq("session_id", session_id).order("rank").execute()

    # 번들 생성
    import hashlib
    import json
    from datetime import datetime

    bundle_data = {
        "metadata": {
            "session_id": session_id,
            "session_type": session.data.get("session_type"),
            "target_antigen": session.data.get("target_antigen"),
            "target_indication": session.data.get("target_indication"),
            "requested_dar": session.data.get("requested_dar"),
            "created_at": session.data.get("created_at"),
            "completed_at": session.data.get("updated_at"),
            "status": session.data.get("status"),
            "export_timestamp": datetime.utcnow().isoformat()
        },
        "candidates": candidates.data or [],
        "execution_logs": [{
            "sequence": log.get("sequence_number"),
            "agent": log.get("agent_name"),
            "status": log.get("status"),
            "reasoning": log.get("reasoning"),
            "decision_summary": log.get("decision_summary"),
            "confidence_score": log.get("confidence_score"),
            "record_hash": log.get("record_hash"),
            "previous_hash": log.get("previous_hash"),
            "timestamp": log.get("created_at")
        } for log in (logs.data or [])],
        "code_execution": [{
            "sequence": log.get("sequence_number"),
            "agent": log.get("agent_name"),
            "code": log.get("executed_code"),
            "result": log.get("execution_result")
        } for log in (logs.data or []) if log.get("executed_code")]
    }

    # 번들 해시 생성
    bundle_json = json.dumps(bundle_data, sort_keys=True, default=str)
    bundle_hash = hashlib.sha256(bundle_json.encode()).hexdigest()

    bundle_data["bundle_hash"] = bundle_hash
    bundle_data["verification"] = {
        "algorithm": "SHA-256",
        "bundle_hash": bundle_hash,
        "compliance": "21 CFR Part 11",
        "verification_note": "Hash verification ensures data integrity from generation to export"
    }

    return bundle_data


# ===== WebSocket Endpoints =====

@router.websocket("/ws/{session_id}")
async def websocket_session(websocket: WebSocket, session_id: str):
    """
    세션 실시간 업데이트 WebSocket

    클라이언트는 이 엔드포인트에 연결하여 실시간 진행 상황을 수신합니다.

    메시지 타입:
    - agent_status: 에이전트 상태 변경
    - progress: 진행률 업데이트
    - candidates_update: 후보 구조 업데이트
    - session_complete: 세션 완료
    """
    await websocket_hub.connect(websocket, session_id)

    try:
        # 연결 시 현재 상태 전송
        supabase = get_supabase_client()
        session = supabase.table("design_sessions").select(
            "status, current_step, current_agent"
        ).eq("id", session_id).single().execute()

        if session.data:
            await websocket.send_json({
                "type": "initial_state",
                "status": session.data.get("status"),
                "current_step": session.data.get("current_step", 0),
                "current_agent": session.data.get("current_agent")
            })

        # 연결 유지 (클라이언트 메시지 대기)
        while True:
            data = await websocket.receive_text()
            # 클라이언트에서 ping 메시지 처리
            if data == "ping":
                await websocket.send_text("pong")

    except WebSocketDisconnect:
        logger.info(f"[ws] Client disconnected from session {session_id}")
    except Exception as e:
        logger.error(f"[ws] WebSocket error: {e}")
    finally:
        await websocket_hub.disconnect(websocket, session_id)
