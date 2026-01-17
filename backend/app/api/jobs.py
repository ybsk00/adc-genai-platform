"""
Jobs API - 시뮬레이션 작업 관리
ADC 분석 작업 생성, 상태 조회, 결과 조회
"""
from fastapi import APIRouter, HTTPException, BackgroundTasks, Depends
from pydantic import BaseModel, Field
from typing import Optional, Dict, Any
from uuid import uuid4
from datetime import datetime

from app.agents import ADCInput, run_adc_analysis, JobStatus, AgentStatus
from app.core.config import settings

router = APIRouter()

# In-memory 작업 저장소 (TODO: Redis/Supabase로 교체)
job_store: Dict[str, Dict[str, Any]] = {}


class SimulationRequest(BaseModel):
    """시뮬레이션 요청 모델"""
    target_name: str = Field(..., description="타겟 이름 (e.g., HER2)")
    antibody_type: str = Field(..., description="항체 타입")
    custom_sequence: Optional[str] = None
    payload_id: str = Field(..., description="약물 ID (e.g., mmae)")
    linker_id: str = Field(..., description="링커 ID (e.g., val-cit)")
    dar: int = Field(default=4, ge=1, le=8)
    mode: str = Field(default="deep", description="fast or deep")
    job_name: Optional[str] = None


class JobResponse(BaseModel):
    """작업 응답 모델"""
    job_id: str
    status: str
    progress: int
    message: Optional[str] = None


class JobStatusResponse(BaseModel):
    """작업 상태 응답 (폴링용)"""
    job_id: str
    status: str
    progress: int
    agents: Dict[str, str]


class JobResultResponse(BaseModel):
    """작업 결과 응답"""
    job_id: str
    status: str
    grade: Optional[str] = None
    recommendation: Optional[str] = None
    summary: Optional[str] = None
    scores: Optional[Dict[str, float]] = None
    structure_analysis: Optional[Dict[str, Any]] = None
    toxicity_risks: Optional[list] = None
    patent_landscape: Optional[list] = None
    competitors: Optional[list] = None
    clinical_protocol: Optional[Dict[str, Any]] = None
    report_url: Optional[str] = None


async def run_analysis_background(job_id: str, input_data: dict, user_id: str):
    """
    백그라운드에서 ADC 분석 실행
    TODO: Celery 태스크로 분리
    """
    try:
        job_store[job_id]["status"] = "running"
        
        # LangGraph 오케스트레이터 실행
        final_state = await run_adc_analysis(input_data, job_id, user_id)
        
        # 결과 저장
        job_store[job_id] = {
            "status": final_state.status.value,
            "progress": 100,
            "grade": final_state.final_grade,
            "recommendation": final_state.recommendation,
            "summary": final_state.executive_summary,
            "scores": final_state.scores,
            "structure_analysis": final_state.structure_analysis.dict() if final_state.structure_analysis else None,
            "toxicity_risks": [r.dict() for r in final_state.toxicity_risks],
            "patent_landscape": [p.dict() for p in final_state.patent_landscape],
            "competitors": [c.dict() for c in final_state.competitors],
            "clinical_protocol": final_state.clinical_protocol,
            "report_url": final_state.report_url,
            "agents": {k: v.status.value for k, v in final_state.agents.items()},
            "errors": final_state.errors,
            "completed_at": datetime.utcnow().isoformat()
        }
        
    except Exception as e:
        job_store[job_id] = {
            "status": "failed",
            "progress": 0,
            "error": str(e)
        }


@router.post("", response_model=JobResponse)
async def create_simulation(
    req: SimulationRequest, 
    background_tasks: BackgroundTasks
):
    """
    시뮬레이션 작업 생성
    
    Returns:
        JobResponse: 생성된 작업 정보
    """
    job_id = f"job_{uuid4().hex[:8]}"
    user_id = "user_demo"  # TODO: JWT에서 추출
    
    # 크레딧 차감 (TODO: 실제 구현)
    credits_required = 10 if req.mode == "deep" else 1
    
    # 초기 상태 저장
    job_store[job_id] = {
        "status": "queued",
        "progress": 0,
        "created_at": datetime.utcnow().isoformat(),
        "input": req.dict(),
        "agents": {
            "structure": "pending",
            "toxicology": "pending",
            "patent": "pending",
            "competitor": "pending",
            "clinical": "pending",
            "report": "pending",
        }
    }
    
    # 백그라운드에서 분석 실행
    input_data = req.dict()
    background_tasks.add_task(run_analysis_background, job_id, input_data, user_id)
    
    return JobResponse(
        job_id=job_id,
        status="queued",
        progress=0,
        message=f"Simulation job submitted. {credits_required} credits will be deducted."
    )


@router.get("/{job_id}/status", response_model=JobStatusResponse)
async def get_job_status(job_id: str):
    """
    작업 상태 조회 (폴링용)
    프론트엔드에서 3초마다 호출
    """
    if job_id not in job_store:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = job_store[job_id]
    
    # 진행률 계산
    progress = job.get("progress", 0)
    if job.get("status") == "running":
        agents = job.get("agents", {})
        done_count = sum(1 for s in agents.values() if s == "done")
        progress = int(done_count / 6 * 100)
    
    return JobStatusResponse(
        job_id=job_id,
        status=job.get("status", "unknown"),
        progress=progress,
        agents=job.get("agents", {})
    )


@router.get("/{job_id}", response_model=JobResultResponse)
async def get_job_result(job_id: str):
    """
    작업 결과 조회
    """
    if job_id not in job_store:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = job_store[job_id]
    
    return JobResultResponse(
        job_id=job_id,
        status=job.get("status", "unknown"),
        grade=job.get("grade"),
        recommendation=job.get("recommendation"),
        summary=job.get("summary"),
        scores=job.get("scores"),
        structure_analysis=job.get("structure_analysis"),
        toxicity_risks=job.get("toxicity_risks"),
        patent_landscape=job.get("patent_landscape"),
        competitors=job.get("competitors"),
        clinical_protocol=job.get("clinical_protocol"),
        report_url=job.get("report_url")
    )


@router.get("/{job_id}/report-url")
async def get_report_download_url(job_id: str):
    """
    리포트 다운로드 URL 생성 (Pre-signed URL, 5분 유효)
    """
    if job_id not in job_store:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = job_store[job_id]
    
    if job.get("status") != "completed":
        raise HTTPException(status_code=400, detail="Job not completed yet")
    
    # TODO: S3 Pre-signed URL 생성
    presigned_url = f"https://s3.amazonaws.com/adc-reports/{job_id}.pdf?X-Amz-Expires=300"
    
    return {
        "url": presigned_url,
        "expires_in": 300
    }
