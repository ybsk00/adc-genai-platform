"""
Design Engine API Routes
ADC 설계 세션 및 워크플로우 API
"""
from fastapi import APIRouter, Depends, HTTPException, BackgroundTasks, WebSocket, WebSocketDisconnect
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from typing import Optional, List, AsyncGenerator
import logging
import asyncio
import json

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


# ===== Auto-Benchmarking API (Phase 2) =====

class AutoBenchmarkRequest(BaseModel):
    """자동 벤치마킹 요청"""
    target_antigen: str
    payload_smiles: str
    linker_type: Optional[str] = None
    linker_smiles: Optional[str] = None
    dar: Optional[float] = None
    our_metrics: Optional[dict] = None


@router.post("/benchmark/auto")
async def run_auto_benchmark(request: AutoBenchmarkRequest):
    """
    자동 벤치마킹 실행

    신규 설계안과 기존 Catalog ADC 3종 자동 비교
    - Tanimoto 유사도 기반 매칭
    - 독성, 결합력, 안정성, 비용 비교
    - 종합 추천 생성
    """
    from app.services.auto_benchmarking_service import get_auto_benchmarking_service

    service = get_auto_benchmarking_service()

    try:
        result = await service.auto_benchmark(
            target_antigen=request.target_antigen,
            payload_smiles=request.payload_smiles,
            linker_type=request.linker_type,
            linker_smiles=request.linker_smiles,
            our_metrics=request.our_metrics,
            dar=request.dar
        )
        return result

    except Exception as e:
        logger.error(f"[auto-benchmark] Error: {e}")
        raise HTTPException(500, f"Auto-benchmarking failed: {str(e)}")


@router.get("/benchmark/session/{session_id}")
async def get_session_benchmark(session_id: str):
    """
    세션의 벤치마킹 결과 조회

    설계 세션에서 생성된 후보의 자동 벤치마킹 결과 반환
    """
    from app.services.auto_benchmarking_service import get_auto_benchmarking_service

    supabase = get_supabase_client()

    # 세션 데이터 조회
    session = supabase.table("design_sessions").select(
        "id, target_antigen, linker_preference, requested_dar"
    ).eq("id", session_id).single().execute()

    if not session.data:
        raise HTTPException(404, "Session not found")

    # 후보 데이터 조회
    candidates = supabase.table("design_candidates").select(
        "smiles, metrics"
    ).eq("session_id", session_id).order("rank").limit(1).execute()

    if not candidates.data:
        raise HTTPException(404, "No candidates found for session")

    top_candidate = candidates.data[0]

    # 자동 벤치마킹 실행
    service = get_auto_benchmarking_service()

    try:
        result = await service.auto_benchmark(
            target_antigen=session.data.get("target_antigen", ""),
            payload_smiles=top_candidate.get("smiles", ""),
            linker_type=session.data.get("linker_preference"),
            our_metrics=top_candidate.get("metrics", {}),
            dar=session.data.get("requested_dar")
        )

        return {
            "session_id": session_id,
            "candidate_smiles": top_candidate.get("smiles", "")[:50] + "...",
            **result
        }

    except Exception as e:
        logger.error(f"[session-benchmark] Error: {e}")
        raise HTTPException(500, f"Session benchmarking failed: {str(e)}")


@router.get("/catalog/search")
async def search_catalog_adcs(
    target: Optional[str] = None,
    linker_type: Optional[str] = None,
    manufacturer: Optional[str] = None,
    limit: int = 20
):
    """
    Catalog ADC 검색

    타겟, 링커 타입, 제조사 기반 검색
    """
    supabase = get_supabase_client()

    query = supabase.table("catalog_adcs").select(
        "id, manufacturer, catalog_number, product_name, "
        "target, target_normalized, linker_type, dar, "
        "price_usd, availability"
    )

    if target:
        query = query.or_(f"target_normalized.ilike.%{target}%,target.ilike.%{target}%")

    if linker_type:
        query = query.eq("linker_type", linker_type)

    if manufacturer:
        query = query.ilike("manufacturer", f"%{manufacturer}%")

    result = query.limit(limit).execute()

    return {
        "total": len(result.data or []),
        "catalog_adcs": result.data or []
    }


# ===== Sandbox Stress Test API (Phase 2) =====

class PhysicalValidationRequest(BaseModel):
    """물리 검증 요청"""
    smiles: str
    molecule_name: Optional[str] = None
    generate_3d: bool = True
    save_to_db: bool = False


@router.post("/validate/physical")
async def validate_physical_structure(request: PhysicalValidationRequest):
    """
    물리적 구조 검증 (Sandbox Stress Test)

    - SMILES 유효성 검증
    - 원자가 규칙 검증
    - 물리화학적 속성 검증
    - 3D 구조 기반 검증 (steric clash, bond length, bond angle)
    - Druglikeness 검증 (ADC 특화)
    """
    from app.services.physical_validator import get_physical_validator

    validator = get_physical_validator()

    try:
        result = await validator.validate_structure(
            smiles=request.smiles,
            molecule_name=request.molecule_name,
            generate_3d=request.generate_3d,
            save_to_db=request.save_to_db
        )
        return result

    except Exception as e:
        logger.error(f"[physical-validation] Error: {e}")
        raise HTTPException(500, f"Physical validation failed: {str(e)}")


@router.post("/validate/session/{session_id}")
async def validate_session_candidates(session_id: str):
    """
    세션의 모든 후보에 대한 물리 검증 실행

    설계된 ADC 후보들의 물리적 타당성 종합 검증
    """
    from app.services.physical_validator import get_physical_validator

    supabase = get_supabase_client()

    # 세션 확인
    session = supabase.table("design_sessions").select("id").eq(
        "id", session_id
    ).single().execute()

    if not session.data:
        raise HTTPException(404, "Session not found")

    # 후보 조회
    candidates = supabase.table("design_candidates").select(
        "id, smiles, rank"
    ).eq("session_id", session_id).order("rank").execute()

    if not candidates.data:
        raise HTTPException(404, "No candidates found for session")

    validator = get_physical_validator()
    results = []

    for candidate in candidates.data:
        try:
            result = await validator.validate_structure(
                smiles=candidate.get("smiles", ""),
                session_id=session_id,
                molecule_name=f"Candidate Rank {candidate.get('rank', 'N/A')}",
                generate_3d=True,
                save_to_db=True
            )
            results.append({
                "candidate_id": candidate["id"],
                "rank": candidate.get("rank"),
                "overall_status": result.get("overall_status"),
                "summary": result.get("summary")
            })
        except Exception as e:
            results.append({
                "candidate_id": candidate["id"],
                "rank": candidate.get("rank"),
                "overall_status": "error",
                "error": str(e)
            })

    # 전체 요약
    total = len(results)
    passed = sum(1 for r in results if r.get("overall_status") == "pass")
    warnings = sum(1 for r in results if r.get("overall_status") == "warning")
    failed = sum(1 for r in results if r.get("overall_status") in ["fail", "error"])

    return {
        "session_id": session_id,
        "total_candidates": total,
        "summary": {
            "passed": passed,
            "with_warnings": warnings,
            "failed": failed
        },
        "results": results
    }


@router.get("/validate/session/{session_id}/history")
async def get_validation_history(session_id: str):
    """
    세션의 물리 검증 이력 조회
    """
    supabase = get_supabase_client()

    result = supabase.table("physical_validations").select(
        "check_name, check_category, passed, severity, actual_value, details, created_at"
    ).eq("session_id", session_id).order("created_at", desc=True).limit(100).execute()

    # 카테고리별 그룹화
    by_category = {}
    for v in result.data or []:
        cat = v.get("check_category", "other")
        if cat not in by_category:
            by_category[cat] = []
        by_category[cat].append(v)

    return {
        "session_id": session_id,
        "total_validations": len(result.data or []),
        "by_category": by_category
    }


# ===== One-Click ADC Navigator API (Phase 3) =====

class NavigatorRequest(BaseModel):
    """Navigator 요청"""
    disease_name: str
    session_id: Optional[str] = None


class DiseaseSuggestion(BaseModel):
    """질환명 자동완성 항목"""
    disease: str
    primary_targets: List[str]
    approved_adc_count: int
    avg_orr: Optional[float] = None


@router.post("/navigator/run")
async def run_navigator(request: NavigatorRequest, background_tasks: BackgroundTasks):
    """
    One-Click ADC Navigator 실행

    질환명 하나만 입력하면 6인 에이전트가 협업하여
    최적의 ADC 설계안을 자동 생성합니다.

    3단계 파이프라인:
    - Step 1: Target & Antibody Match (Librarian)
    - Step 2: Linker & Payload Coupling (Alchemist)
    - Step 3: Simulation & Audit (Coder + Auditor)
    """
    from app.agents.navigator_orchestrator import get_navigator_orchestrator
    import uuid

    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    # 세션 ID 생성 또는 사용
    session_id = request.session_id or str(uuid.uuid4())

    try:
        # 세션 생성
        supabase.table("navigator_sessions").insert({
            "id": session_id,
            "user_id": user_id,
            "disease_name": request.disease_name,
            "status": "running",
            "current_step": 0,
            "total_steps": 5
        }).execute()

        # 백그라운드에서 Navigator 실행
        async def run_navigator_background():
            from app.agents.navigator_orchestrator import run_one_click_navigator
            try:
                result = await run_one_click_navigator(
                    disease_name=request.disease_name,
                    session_id=session_id,
                    user_id=user_id
                )

                # 결과 저장
                supabase.table("navigator_sessions").update({
                    "status": "completed",
                    "antibody_candidates": [
                        {
                            "id": ab.antibody_id,
                            "name": ab.name,
                            "target_protein": ab.target_protein,
                            "clinical_score": ab.clinical_score,
                            "match_confidence": ab.match_confidence
                        }
                        for ab in result.antibody_candidates
                    ],
                    "primary_target": result.antibody_candidates[0].target_protein if result.antibody_candidates else None,
                    "golden_combination": {
                        "antibody": result.golden_combination.antibody.name,
                        "linker": {
                            "type": result.golden_combination.linker.type,
                            "smiles": result.golden_combination.linker.smiles
                        },
                        "payload": {
                            "class": result.golden_combination.payload.class_name,
                            "smiles": result.golden_combination.payload.smiles
                        },
                        "dar": result.golden_combination.dar,
                        "historical_orr": result.golden_combination.historical_performance.get("orr_pct") if result.golden_combination.historical_performance else None
                    },
                    "physics_verified": result.physics_verified,
                    "virtual_trial": {
                        "predicted_orr": result.virtual_trial.predicted_orr,
                        "predicted_pfs_months": result.virtual_trial.predicted_pfs_months,
                        "predicted_os_months": result.virtual_trial.predicted_os_months,
                        "pk_data": result.virtual_trial.pk_data,
                        "tumor_data": result.virtual_trial.tumor_data,
                        "confidence": result.virtual_trial.confidence
                    },
                    "predicted_orr": result.virtual_trial.predicted_orr,
                    "lineage_data": result.digital_lineage,
                    "current_step": 5
                }).eq("id", session_id).execute()

            except Exception as e:
                logger.exception(f"[navigator-bg] Error: {e}")
                supabase.table("navigator_sessions").update({
                    "status": "failed",
                    "error_message": str(e)
                }).eq("id", session_id).execute()

        # 백그라운드 태스크 추가
        background_tasks.add_task(run_navigator_background)

        return {
            "success": True,
            "session_id": session_id,
            "status": "running",
            "message": "Navigator started. Connect to /api/design/navigator/stream/{session_id} for progress."
        }

    except Exception as e:
        logger.exception(f"[navigator] Error: {e}")
        raise HTTPException(500, f"Navigator failed to start: {str(e)}")


@router.get("/navigator/session/{session_id}")
async def get_navigator_session(session_id: str):
    """
    Navigator 세션 상태 조회
    """
    supabase = get_supabase_client()

    result = supabase.table("navigator_sessions").select("*").eq(
        "id", session_id
    ).single().execute()

    if not result.data:
        raise HTTPException(404, "Navigator session not found")

    return result.data


@router.get("/navigator/sessions")
async def list_navigator_sessions(
    status: Optional[str] = None,
    limit: int = 20
):
    """
    Navigator 세션 목록 조회
    """
    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    query = supabase.table("navigator_sessions").select(
        "id, disease_name, status, physics_verified, predicted_orr, "
        "current_step, total_steps, created_at, completed_at"
    ).eq("user_id", user_id)

    if status:
        query = query.eq("status", status)

    result = query.order("created_at", desc=True).limit(limit).execute()

    return {
        "total": len(result.data or []),
        "sessions": result.data or []
    }


@router.get("/navigator/suggest-disease")
async def suggest_diseases(query: str, limit: int = 10):
    """
    질환명 자동완성

    질환명의 일부를 입력하면 관련 질환과 주요 타겟을 반환합니다.
    """
    if len(query) < 2:
        return {"suggestions": []}

    supabase = get_supabase_client()

    try:
        result = supabase.rpc("suggest_diseases", {
            "query_text": query,
            "max_results": limit
        }).execute()

        suggestions = []
        for row in result.data or []:
            suggestions.append({
                "disease": row.get("disease"),
                "primary_targets": row.get("primary_targets", []),
                "approved_adc_count": row.get("approved_adc_count", 0),
                "avg_orr": row.get("avg_orr")
            })

        return {"suggestions": suggestions}

    except Exception as e:
        logger.warning(f"[suggest-disease] Error: {e}")
        # Fallback: 직접 검색
        result = supabase.table("antibody_library").select(
            "related_disease, target_normalized"
        ).ilike("related_disease", f"%{query}%").limit(limit).execute()

        disease_map = {}
        for row in result.data or []:
            disease = row.get("related_disease")
            if disease:
                if disease not in disease_map:
                    disease_map[disease] = []
                target = row.get("target_normalized")
                if target and target not in disease_map[disease]:
                    disease_map[disease].append(target)

        suggestions = [
            {
                "disease": disease,
                "primary_targets": targets[:3],
                "approved_adc_count": 0,
                "avg_orr": None
            }
            for disease, targets in disease_map.items()
        ]

        return {"suggestions": suggestions}


@router.get("/navigator/disease-targets/{disease}")
async def get_disease_targets(disease: str):
    """
    특정 질환에 대한 타겟 단백질 목록 조회
    """
    supabase = get_supabase_client()

    # antibody_library에서 검색
    ab_result = supabase.table("antibody_library").select(
        "target_protein, target_normalized"
    ).ilike("related_disease", f"%{disease}%").execute()

    # golden_set에서 검색
    gs_result = supabase.table("golden_set").select(
        "target_1, target_2, clinical_status, orr_pct"
    ).ilike("indication", f"%{disease}%").execute()

    # 타겟 집계
    target_stats = {}
    for row in ab_result.data or []:
        target = row.get("target_normalized") or row.get("target_protein")
        if target:
            if target not in target_stats:
                target_stats[target] = {"count": 0, "max_orr": None, "approved": False}
            target_stats[target]["count"] += 1

    for row in gs_result.data or []:
        for target_key in ["target_1", "target_2"]:
            target = row.get(target_key)
            if target:
                if target not in target_stats:
                    target_stats[target] = {"count": 0, "max_orr": None, "approved": False}
                orr = row.get("orr_pct")
                if orr and (target_stats[target]["max_orr"] is None or orr > target_stats[target]["max_orr"]):
                    target_stats[target]["max_orr"] = orr
                if row.get("clinical_status") == "Approved":
                    target_stats[target]["approved"] = True

    # 정렬 (approved > max_orr > count)
    sorted_targets = sorted(
        target_stats.items(),
        key=lambda x: (x[1]["approved"], x[1]["max_orr"] or 0, x[1]["count"]),
        reverse=True
    )

    return {
        "disease": disease,
        "targets": [
            {
                "target": target,
                "antibody_count": stats["count"],
                "max_orr": stats["max_orr"],
                "has_approved_adc": stats["approved"]
            }
            for target, stats in sorted_targets[:10]
        ]
    }


@router.get("/navigator/stream/{session_id}")
async def stream_navigator_progress(session_id: str):
    """
    Navigator 진행 상황 SSE 스트리밍

    Server-Sent Events를 통해 실시간 진행 상황 전송
    """
    async def event_generator() -> AsyncGenerator[str, None]:
        supabase = get_supabase_client()
        last_step = 0
        max_wait = 120  # 최대 2분 대기
        elapsed = 0

        while elapsed < max_wait:
            try:
                # 세션 상태 조회
                result = supabase.table("navigator_sessions").select(
                    "status, current_step, error_message"
                ).eq("id", session_id).single().execute()

                if not result.data:
                    yield f"data: {json.dumps({'type': 'error', 'message': 'Session not found'})}\n\n"
                    break

                session = result.data
                current_step = session.get("current_step", 0)
                status = session.get("status", "pending")

                # 단계 변경 시 이벤트 전송
                if current_step > last_step:
                    step_names = [
                        "Target & Antibody Match",
                        "Golden Combination",
                        "Property Calculation",
                        "Physical Validation",
                        "Virtual Trial"
                    ]

                    # 이전 단계 완료 이벤트
                    if last_step > 0:
                        yield f"data: {json.dumps({'type': 'step_complete', 'step': last_step, 'message': f'{step_names[last_step-1]} completed'})}\n\n"

                    # 새 단계 시작 이벤트
                    if current_step <= 5:
                        yield f"data: {json.dumps({'type': 'step_start', 'step': current_step, 'message': f'Running {step_names[current_step-1]}...'})}\n\n"

                    last_step = current_step

                # 완료 상태 확인
                if status == "completed":
                    # 마지막 단계 완료
                    if last_step < 5:
                        for step in range(last_step + 1, 6):
                            yield f"data: {json.dumps({'type': 'step_complete', 'step': step, 'message': 'Completed'})}\n\n"

                    # 전체 결과 조회
                    full_result = supabase.table("navigator_sessions").select("*").eq(
                        "id", session_id
                    ).single().execute()

                    if full_result.data:
                        # 결과 데이터 정리
                        data = full_result.data
                        result_payload = {
                            "type": "complete",
                            "result": {
                                "session_id": session_id,
                                "disease_name": data.get("disease_name"),
                                "target_protein": data.get("primary_target"),
                                "antibody_candidates": data.get("antibody_candidates", []),
                                "golden_combination": data.get("golden_combination", {}),
                                "physics_verified": data.get("physics_verified", False),
                                "virtual_trial": data.get("virtual_trial", {}),
                                "execution_time_seconds": 0
                            }
                        }
                        yield f"data: {json.dumps(result_payload)}\n\n"
                    break

                elif status == "failed":
                    yield f"data: {json.dumps({'type': 'error', 'message': session.get('error_message', 'Unknown error')})}\n\n"
                    break

                # 1초 대기
                await asyncio.sleep(1)
                elapsed += 1

            except Exception as e:
                logger.error(f"[navigator-stream] Error: {e}")
                yield f"data: {json.dumps({'type': 'error', 'message': str(e)})}\n\n"
                break

        # 타임아웃
        if elapsed >= max_wait:
            yield f"data: {json.dumps({'type': 'error', 'message': 'Timeout waiting for navigator completion'})}\n\n"

    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no"
        }
    )
