"""
Design Runs API - Design Run 관리 엔드포인트
Eng-Fit v0.2: frozen_params 스냅샷으로 재현성 보장
"""
from fastapi import APIRouter, HTTPException, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime

from app.core.supabase import supabase
from app.services.design_run_service import design_run_service

router = APIRouter()


# ============================================================
# Request/Response Models
# ============================================================

class CreateDesignRunRequest(BaseModel):
    """Design Run 생성 요청"""
    project_id: str = Field(..., description="프로젝트 ID")


class DesignRunResponse(BaseModel):
    """Design Run 응답"""
    id: str
    project_id: str
    status: str
    frozen_params: Dict[str, Any]
    created_at: str
    updated_at: str


class UpdateStatusRequest(BaseModel):
    """상태 업데이트 요청"""
    status: str = Field(..., pattern="^(draft|computing|completed|failed)$")


# ============================================================
# API Endpoints
# ============================================================

@router.post("", response_model=DesignRunResponse)
async def create_design_run(req: CreateDesignRunRequest):
    """
    새 Design Run 생성
    - frozen_params에 현재 스코링 파라미터 스냅샷 저장
    - 재현성 및 감사(Audit) 가능
    """
    try:
        result = await design_run_service.create_design_run(req.project_id)
        return DesignRunResponse(**result)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("")
async def list_design_runs(project_id: Optional[str] = None, limit: int = 50):
    """
    Design Run 목록 조회
    - project_id로 필터링 가능
    """
    try:
        query = supabase.table("design_runs").select("*").order("created_at", desc=True).limit(limit)
        
        if project_id:
            query = query.eq("project_id", project_id)
        
        res = query.execute()
        return res.data or []
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{run_id}")
async def get_design_run(run_id: str):
    """
    특정 Design Run 상세 조회
    - frozen_params 포함
    """
    try:
        result = await design_run_service.get_run(run_id)
        if not result:
            raise HTTPException(status_code=404, detail="Design Run not found")
        return result
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.patch("/{run_id}/status")
async def update_design_run_status(run_id: str, req: UpdateStatusRequest):
    """
    Design Run 상태 업데이트
    """
    try:
        await design_run_service.update_run_status(run_id, req.status)
        return {"status": "updated", "run_id": run_id, "new_status": req.status}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{run_id}/frozen-params")
async def get_frozen_params(run_id: str):
    """
    Design Run의 frozen_params 조회 (감사용)
    """
    try:
        result = await design_run_service.get_run(run_id)
        if not result:
            raise HTTPException(status_code=404, detail="Design Run not found")
        return {
            "run_id": run_id,
            "frozen_params": result.get("frozen_params", {}),
            "frozen_at": result.get("frozen_params", {}).get("frozen_at")
        }
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
