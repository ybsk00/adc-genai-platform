"""
Admin API - 관리자 전용 엔드포인트
KPI 조회, 크레딧 지급, Golden Set 관리 (Human-in-the-Loop), 프롬프트 관리
PDF 업로드, PubMed 크롤러 트리거
"""
from fastapi import APIRouter, HTTPException, Depends, UploadFile, File, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Optional, List
from datetime import datetime
from uuid import uuid4
import tempfile
import os

from app.core.config import settings
from app.core.supabase import supabase
from app.services.pdf_parser import pdf_parser_service
from app.services.pubmed_crawler import pubmed_crawler

router = APIRouter()


# ============================================================
# Request/Response Models
# ============================================================

class CreditGrantRequest(BaseModel):
    """크레딧 지급 요청"""
    user_id: str
    amount: int = Field(..., gt=0)
    reason: str = Field(..., min_length=5, description="지급 사유 (매출 정산용)")


class CreditGrantResponse(BaseModel):
    """크레딧 지급 응답"""
    status: str
    user_id: str
    amount: int
    new_balance: int
    transaction_id: str


class GoldenSetDraft(BaseModel):
    """Golden Set 검토 대기 항목"""
    id: str
    drug_name: str
    target: str
    payload: str
    linker: Optional[str]
    enrichment_source: str
    created_at: str


class ApproveRequest(BaseModel):
    """승인 요청"""
    note: Optional[str] = None


class RejectRequest(BaseModel):
    """반려 요청"""
    reason: str = Field(..., min_length=5, description="반려 사유")


class PromptUpdateRequest(BaseModel):
    """프롬프트 업데이트 요청"""
    system_prompt: str
    version: Optional[str] = None  # 없으면 자동 생성


# ============================================================
# KPI Dashboard
# ============================================================

@router.get("/stats")
async def get_admin_stats():
    """
    대시보드 KPI 조회
    TODO: Supabase에서 실제 데이터 조회
    """
    # TODO: 실제 DB 쿼리
    return {
        "mrr": 12450,
        "total_users": 1204,
        "new_users_today": 15,
        "active_simulations": 45,
        "error_rate": 0.5,
        "pending_drafts": 7,  # 검토 대기 중인 Golden Set
    }


# ============================================================
# Credit Management
# ============================================================

@router.post("/credits/grant", response_model=CreditGrantResponse)
async def grant_credits(req: CreditGrantRequest):
    """
    특정 유저 크레딧 지급
    [Dev Note] 이 기록은 credit_transactions 테이블에 저장되어 매출 정산에 사용됨
    """
    # TODO: Supabase 트랜잭션으로 구현
    # 1. profiles 테이블에서 현재 크레딧 조회
    # 2. 크레딧 업데이트
    # 3. credit_transactions에 기록
    
    transaction_id = f"txn_{uuid4().hex[:8]}"
    
    # Mock response
    return CreditGrantResponse(
        status="success",
        user_id=req.user_id,
        amount=req.amount,
        new_balance=500 + req.amount,  # Mock
        transaction_id=transaction_id
    )


# ============================================================
# Golden Set Management (Human-in-the-Loop)
# ============================================================

@router.get("/goldenset/drafts", response_model=List[GoldenSetDraft])
async def get_golden_set_drafts():
    """
    검토 대기 중인 Golden Set 목록 조회
    
    Worker가 수집한 데이터가 'draft' 상태로 저장되며,
    관리자가 검토 후 승인해야 RAG에 반영됨
    """
    try:
        # Supabase 쿼리: status가 'draft'인 항목 조회, 최신순 정렬
        response = supabase.table("golden_set") \
            .select("*") \
            .eq("status", "draft") \
            .order("created_at", desc=True) \
            .execute()
            
        drafts = []
        for item in response.data:
            drafts.append(GoldenSetDraft(
                id=item.get("id"),
                drug_name=item.get("drug_name") or "Unknown",
                target=item.get("target") or "Unknown",
                payload=item.get("payload"),
                linker=item.get("linker"),
                enrichment_source=item.get("enrichment_source") or "unknown",
                created_at=item.get("created_at")
            ))
            
        return drafts
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/goldenset/{golden_set_id}")
async def get_golden_set_detail(golden_set_id: str):
    """
    Golden Set 상세 조회 (raw_data 포함)
    관리자가 검토 시 원본 데이터 확인용
    """
    try:
        response = supabase.table("golden_set") \
            .select("*") \
            .eq("id", golden_set_id) \
            .execute()
            
        if not response.data:
            raise HTTPException(status_code=404, detail="Golden Set not found")
            
        return response.data[0]
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/goldenset/{golden_set_id}/approve")
async def approve_golden_set(golden_set_id: str, req: ApproveRequest):
    """
    Golden Set 승인 및 RAG 인덱싱
    
    Logic:
    1. DB의 status를 'approved'로 변경
    2. (TODO) 해당 데이터의 description 텍스트를 OpenAI Embedding API로 벡터화 및 저장
    """
    try:
        # 1. 상태 업데이트
        update_data = {
            "status": "approved",
            "approved_at": datetime.utcnow().isoformat()
            # reviewer_id는 현재 인증된 사용자 ID를 넣어야 하지만, 일단 생략
        }
        
        response = supabase.table("golden_set") \
            .update(update_data) \
            .eq("id", golden_set_id) \
            .execute()
            
        if not response.data:
            raise HTTPException(status_code=404, detail="Golden Set not found")
            
        # TODO: RAG 인덱싱 로직 추가 (Phase 5.2)
        # 현재는 상태 변경만 수행
    
        return {
            "status": "approved",
            "golden_set_id": golden_set_id,
            "message": "Data approved (RAG indexing pending)",
            "indexed_at": datetime.utcnow().isoformat()
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/goldenset/{golden_set_id}/reject")
async def reject_golden_set(golden_set_id: str, req: RejectRequest):
    """
    Golden Set 반려
    
    Logic:
    - status를 'rejected'로 변경
    - reviewer_note에 반려 사유 기록
    """
    try:
        update_data = {
            "status": "rejected",
            "reviewer_note": req.reason,
            "rejected_at": datetime.utcnow().isoformat()
        }
        
        response = supabase.table("golden_set") \
            .update(update_data) \
            .eq("id", golden_set_id) \
            .execute()
            
        if not response.data:
            raise HTTPException(status_code=404, detail="Golden Set not found")
    
        return {
            "status": "rejected",
            "golden_set_id": golden_set_id,
            "reason": req.reason,
            "message": "Data rejected"
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/goldenset/sync")
async def trigger_golden_set_sync():
    """
    외부 데이터 소스 동기화 트리거
    Worker가 PubMed/ClinicalTrials 데이터를 수집하고 draft로 저장
    """
    # TODO: Celery 태스크 또는 Cloud Run Job 트리거
    
    return {
        "status": "started",
        "message": "Data sync job triggered. New data will appear in drafts.",
        "job_id": f"sync_{uuid4().hex[:8]}"
    }


@router.post("/goldenset/reindex")
async def reindex_rag():
    """
    전체 RAG 재인덱싱
    모든 approved 데이터의 임베딩을 재생성
    
    [WARNING] 비용과 시간이 많이 소요됨
    """
    # TODO: 백그라운드 태스크로 실행
    
    return {
        "status": "started",
        "message": "Re-indexing all approved golden set data",
        "estimated_time": "30 minutes",
        "estimated_cost": "$15"
    }


# ============================================================
# Prompt Management
# ============================================================

@router.get("/prompts")
async def get_all_prompts():
    """모든 에이전트의 현재 프롬프트 조회"""
    # TODO: Supabase 쿼리
    
    return {
        "structure": {"version": "v1.2", "is_live": True},
        "toxicology": {"version": "v1.0", "is_live": True, "draft": "v1.1"},
        "patent": {"version": "v1.0", "is_live": True},
        "competitor": {"version": "v1.0", "is_live": True},
        "clinical": {"version": "v1.0", "is_live": True},
        "report": {"version": "v2.0", "is_live": True},
    }


@router.get("/prompts/{agent_id}")
async def get_prompt(agent_id: str):
    """특정 에이전트의 프롬프트 조회 (버전 히스토리 포함)"""
    # TODO: Supabase 쿼리
    
    return {
        "agent_id": agent_id,
        "current_version": "v1.0",
        "system_prompt": "You are a Structure Analysis Agent...",
        "versions": [
            {"version": "v1.0", "is_live": True, "created_at": "2026-01-10"},
            {"version": "v0.9", "is_live": False, "created_at": "2026-01-05"},
        ]
    }


@router.put("/prompts/{agent_id}")
async def update_prompt(agent_id: str, req: PromptUpdateRequest):
    """
    에이전트 프롬프트 수정 (Draft 저장)
    Live 배포는 별도 publish 엔드포인트에서 처리
    """
    # TODO: agent_prompts 테이블에 새 버전 INSERT
    
    version = req.version or f"v{datetime.now().strftime('%Y%m%d%H%M')}"
    
    return {
        "agent_id": agent_id,
        "version": version,
        "status": "draft_saved",
        "message": f"Prompt draft {version} saved. Use /publish to go live."
    }


@router.post("/prompts/{agent_id}/publish")
async def publish_prompt(agent_id: str, version: str):
    """
    프롬프트 Live 배포
    기존 Live 버전은 is_live=False로 변경
    """
    # TODO: Transaction으로 처리
    # 1. UPDATE agent_prompts SET is_live = FALSE WHERE agent_id = ? AND is_live = TRUE
    # 2. UPDATE agent_prompts SET is_live = TRUE, published_at = NOW() WHERE agent_id = ? AND version = ?
    
    return {
        "agent_id": agent_id,
        "version": version,
        "status": "published",
        "message": f"Prompt {version} is now live"
    }


# ============================================================
# PDF Upload & RAG Ingestion (Phase 3)
# ============================================================

class UploadResponse(BaseModel):
    """PDF 업로드 응답"""
    status: str
    document_id: Optional[str] = None
    parsed_data: Optional[dict] = None
    error: Optional[str] = None


@router.post("/upload", response_model=UploadResponse)
async def upload_document(
    file: UploadFile = File(...),
    background_tasks: BackgroundTasks = None
):
    """
    PDF 문서 업로드 및 Gemini 멀티모달 파싱
    
    지원 형식: PDF
    Gemini Vision으로 테이블/그래프 포함 추출
    """
    if not file.filename.lower().endswith('.pdf'):
        raise HTTPException(status_code=400, detail="Only PDF files are supported")
    
    try:
        # 임시 파일 저장
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tmp:
            content = await file.read()
            tmp.write(content)
            tmp_path = tmp.name
        
        # Gemini로 PDF 파싱
        result = await pdf_parser_service.parse_pdf(tmp_path)
        
        # 임시 파일 삭제
        os.unlink(tmp_path)
        
        if result.get("status") == "success":
            doc_id = f"doc_{uuid4().hex[:8]}"
            data = result.get("data", {})
            
            # 1. DOI 중복 체크 (Phase 3.3)
            doi = data.get("doi")
            if doi:
                existing = supabase.table("golden_set").select("id").eq("raw_data->>doi", doi).execute()
                if existing.data:
                    return UploadResponse(
                        status="skipped",
                        error="Duplicate DOI found",
                        document_id=existing.data[0]['id']
                    )
            
            # 2. DB 저장 (Draft 상태)
            drug_info = data.get("drug_info", {})
            
            new_record = {
                "drug_name": drug_info.get("drug_name") or "Unknown",
                "target": drug_info.get("target") or "Unknown",
                "antibody": drug_info.get("antibody"),
                "payload": drug_info.get("payload"),
                "linker": drug_info.get("linker"),
                "dar": float(drug_info.get("dar")) if drug_info.get("dar") and str(drug_info.get("dar")).replace('.','',1).isdigit() else None,
                "status": "draft",
                "enrichment_source": "pdf_upload",
                "raw_data": data
            }
            
            # Supabase Insert
            db_res = supabase.table("golden_set").insert(new_record).execute()
            
            return UploadResponse(
                status="success",
                document_id=db_res.data[0]['id'] if db_res.data else doc_id,
                parsed_data=data
            )
        else:
            return UploadResponse(
                status="error",
                error=result.get("error")
            )
            
    except Exception as e:
        return UploadResponse(
            status="error",
            error=str(e)
        )


# ============================================================
# PubMed Crawler (Phase 4)
# ============================================================

class CrawlRequest(BaseModel):
    """크롤러 요청"""
    query: Optional[str] = None
    max_results: int = Field(default=20, le=100)
    days_back: int = Field(default=7, le=30)


@router.post("/crawler/trigger")
async def trigger_pubmed_crawler(req: CrawlRequest):
    """
    PubMed 크롤러 트리거
    
    ADC 관련 논문을 수집하고 AI Gatekeeper로 필터링
    필터링 통과한 논문은 draft 상태로 저장
    """
    try:
        result = await pubmed_crawler.crawl_and_filter(
            query=req.query,
            max_results=req.max_results,
            days_back=req.days_back
        )
        
        saved_count = 0
        skipped_count = 0
        
        # 필터링 통과한 논문을 DB에 draft로 저장
        for article in result.get("relevant", []):
            doi = article.get("doi")
            
            # 1. 중복 체크 (Phase 4.3)
            if doi:
                existing = supabase.table("golden_set").select("id").eq("raw_data->>doi", doi).execute()
                if existing.data:
                    skipped_count += 1
                    continue
            
            # 2. DB 저장
            new_record = {
                "drug_name": "Unknown (Auto-crawled)",
                "target": "Unknown",
                "antibody": None,
                "payload": None,
                "status": "draft",
                "enrichment_source": "pubmed_crawler",
                "raw_data": article
            }
            
            supabase.table("golden_set").insert(new_record).execute()
            saved_count += 1
        
        return {
            "status": "success",
            "job_id": f"crawl_{uuid4().hex[:8]}",
            "total_found": result.get("total_found"),
            "processed": result.get("total_processed"),
            "saved": saved_count,
            "skipped_duplicate": skipped_count,
            "details": result
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/crawler/status/{job_id}")
async def get_crawler_status(job_id: str):
    """크롤러 작업 상태 조회"""
    # TODO: Celery/Redis에서 작업 상태 조회
    return {
        "job_id": job_id,
        "status": "completed",
        "message": "Crawler job status (mock)"
    }


# ============================================================
# Health Check
# ============================================================

@router.get("/health")
async def admin_health():
    """관리자 API 헬스 체크"""
    return {"status": "healthy", "service": "admin-api"}

