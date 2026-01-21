"""
Assay Results API - Assay 결과 관리 엔드포인트
Closed-loop: acceptance_criteria와 함께 유연한 성공 판정
"""
from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any

from app.core.supabase import supabase
from app.services.assay_result_service import assay_result_service

router = APIRouter()


# ============================================================
# Request/Response Models
# ============================================================

class CreateAssayResultRequest(BaseModel):
    """Assay 결과 생성 요청"""
    run_id: str = Field(..., description="Design Run ID")
    molecule_id: str = Field(..., description="분자/후보물질 ID")
    assay_type: str = Field(..., description="분석 유형 (aggregation, binding, cytotoxicity)")
    raw_data: Dict[str, Any] = Field(..., description="실험 데이터")
    acceptance_criteria: Dict[str, Any] = Field(..., description="판정 기준")


class AssayResultResponse(BaseModel):
    """Assay 결과 응답"""
    id: str
    run_id: str
    molecule_id: str
    assay_type: str
    raw_data: Dict[str, Any]
    is_success: bool
    acceptance_criteria: Dict[str, Any]
    confidence_score: float
    created_at: str


class AssaySummary(BaseModel):
    """Assay 결과 요약"""
    total: int
    success_count: int
    failure_count: int
    avg_confidence: float


# ============================================================
# API Endpoints
# ============================================================

@router.post("", response_model=AssayResultResponse)
async def create_assay_result(req: CreateAssayResultRequest):
    """
    Assay 결과 저장
    - acceptance_criteria와 함께 결과 저장
    - is_success는 기준에 따라 자동 판정
    - confidence_score 자동 계산
    """
    try:
        result = await assay_result_service.save_assay_result(
            run_id=req.run_id,
            molecule_id=req.molecule_id,
            assay_type=req.assay_type,
            raw_data=req.raw_data,
            acceptance_criteria=req.acceptance_criteria
        )
        return AssayResultResponse(**result)
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/by-run/{run_id}")
async def get_results_by_run(run_id: str):
    """
    특정 Design Run의 모든 Assay 결과 조회
    """
    try:
        results = await assay_result_service.get_results_by_run(run_id)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/by-molecule/{molecule_id}")
async def get_results_by_molecule(molecule_id: str):
    """
    특정 분자의 모든 Assay 결과 조회
    """
    try:
        results = await assay_result_service.get_results_by_molecule(molecule_id)
        return results
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/summary/{run_id}")
async def get_assay_summary(run_id: str):
    """
    Design Run의 Assay 결과 요약
    """
    try:
        results = await assay_result_service.get_results_by_run(run_id)
        
        if not results:
            return AssaySummary(total=0, success_count=0, failure_count=0, avg_confidence=0.0)
        
        total = len(results)
        success_count = sum(1 for r in results if r.get("is_success"))
        failure_count = total - success_count
        
        confidences = [r.get("confidence_score", 0.5) for r in results]
        avg_confidence = sum(confidences) / len(confidences) if confidences else 0.0
        
        return AssaySummary(
            total=total,
            success_count=success_count,
            failure_count=failure_count,
            avg_confidence=round(avg_confidence, 3)
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/{result_id}")
async def get_assay_result(result_id: str):
    """
    특정 Assay 결과 상세 조회
    """
    try:
        res = supabase.table("assay_results").select("*").eq("id", result_id).execute()
        if not res.data:
            raise HTTPException(status_code=404, detail="Assay Result not found")
        return res.data[0]
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
