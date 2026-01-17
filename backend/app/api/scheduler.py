"""
Data Worker API - 외부 데이터 소스 동기화
ClinicalTrials.gov, PubMed 등에서 데이터 수집 후 draft 상태로 저장
"""
from fastapi import APIRouter, BackgroundTasks, HTTPException
from pydantic import BaseModel
from typing import Optional, List
from datetime import datetime
from uuid import uuid4

from app.services.perplexity_service import enrich_drug_data
from app.core.config import settings

router = APIRouter()


# In-memory 동기화 상태 추적 (TODO: Redis로 대체)
sync_jobs: dict = {}


class SyncJobResponse(BaseModel):
    """동기화 작업 응답"""
    job_id: str
    status: str
    message: str


class SyncJobStatus(BaseModel):
    """동기화 작업 상태"""
    job_id: str
    status: str
    source: str
    records_found: int
    records_drafted: int
    started_at: str
    completed_at: Optional[str] = None
    errors: List[str] = []


async def process_clinical_trials_data(job_id: str):
    """
    ClinicalTrials.gov 데이터 수집 및 처리
    
    [Human-in-the-Loop] 모든 데이터는 'draft' 상태로 저장
    """
    sync_jobs[job_id]["status"] = "running"
    
    try:
        # TODO: ClinicalTrials.gov API 호출
        # Mock data
        raw_records = [
            {
                "nct_id": "NCT12345678",
                "drug_name": "ADC-LIV1-001",
                "target": "LIV-1",
                "antibody": "Anti-LIV-1 mAb",
                "payload": None,  # 비어있음 -> Perplexity로 보강
                "linker": None,   # 비어있음 -> Perplexity로 보강
                "phase": "Phase 1"
            },
            {
                "nct_id": "NCT87654321",
                "drug_name": "HER2-DXd-Candidate",
                "target": "HER2",
                "antibody": "Trastuzumab",
                "payload": "DXd",
                "linker": "GGFG",
                "phase": "Phase 2"
            }
        ]
        
        sync_jobs[job_id]["records_found"] = len(raw_records)
        drafted = 0
        
        for record in raw_records:
            try:
                # Perplexity로 데이터 보강
                enriched = await enrich_drug_data(record)
                
                # TODO: Supabase에 draft 상태로 저장
                # supabase.table("golden_set").insert({
                #     **enriched,
                #     "status": "draft",  # [IMPORTANT] draft로 저장!
                #     "raw_data": record
                # }).execute()
                
                drafted += 1
                sync_jobs[job_id]["records_drafted"] = drafted
                
            except Exception as e:
                sync_jobs[job_id]["errors"].append(f"{record.get('drug_name')}: {str(e)}")
        
        sync_jobs[job_id]["status"] = "completed"
        sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        
        # TODO: data_sync_logs 테이블에 기록
        
    except Exception as e:
        sync_jobs[job_id]["status"] = "failed"
        sync_jobs[job_id]["errors"].append(str(e))


async def process_pubmed_data(job_id: str):
    """
    PubMed/BioRxiv 논문 데이터 수집
    """
    sync_jobs[job_id]["status"] = "running"
    
    try:
        # TODO: PubMed API 호출
        # Mock implementation
        
        sync_jobs[job_id]["records_found"] = 150
        sync_jobs[job_id]["records_drafted"] = 142
        sync_jobs[job_id]["status"] = "completed"
        sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        
    except Exception as e:
        sync_jobs[job_id]["status"] = "failed"
        sync_jobs[job_id]["errors"].append(str(e))


@router.post("/sync/clinical", response_model=SyncJobResponse)
async def sync_clinical_trials(background_tasks: BackgroundTasks):
    """
    ClinicalTrials.gov 데이터 동기화 트리거
    백그라운드에서 실행되며, 수집된 데이터는 draft 상태로 저장
    """
    job_id = f"sync_clinical_{uuid4().hex[:8]}"
    
    sync_jobs[job_id] = {
        "status": "queued",
        "source": "clinical_trials",
        "records_found": 0,
        "records_drafted": 0,
        "started_at": datetime.utcnow().isoformat(),
        "completed_at": None,
        "errors": []
    }
    
    background_tasks.add_task(process_clinical_trials_data, job_id)
    
    return SyncJobResponse(
        job_id=job_id,
        status="queued",
        message="ClinicalTrials.gov sync started. Data will be saved as draft."
    )


@router.post("/sync/pubmed", response_model=SyncJobResponse)
async def sync_pubmed(background_tasks: BackgroundTasks):
    """PubMed/BioRxiv 데이터 동기화 트리거"""
    job_id = f"sync_pubmed_{uuid4().hex[:8]}"
    
    sync_jobs[job_id] = {
        "status": "queued",
        "source": "pubmed",
        "records_found": 0,
        "records_drafted": 0,
        "started_at": datetime.utcnow().isoformat(),
        "completed_at": None,
        "errors": []
    }
    
    background_tasks.add_task(process_pubmed_data, job_id)
    
    return SyncJobResponse(
        job_id=job_id,
        status="queued",
        message="PubMed sync started. Data will be saved as draft."
    )


@router.get("/sync/{job_id}", response_model=SyncJobStatus)
async def get_sync_status(job_id: str):
    """동기화 작업 상태 조회"""
    if job_id not in sync_jobs:
        raise HTTPException(status_code=404, detail="Sync job not found")
    
    job = sync_jobs[job_id]
    
    return SyncJobStatus(
        job_id=job_id,
        **job
    )


@router.get("/health")
async def scheduler_health():
    """스케줄러 API 헬스 체크"""
    return {
        "status": "healthy",
        "active_jobs": len([j for j in sync_jobs.values() if j["status"] == "running"])
    }
