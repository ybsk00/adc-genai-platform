"""
Admin API - 관리자 전용 엔드포인트
KPI 조회, 크레딧 지급, Golden Set 관리 (Human-in-the-Loop), 프롬프트 관리
PDF 업로드, PubMed 크롤러 트리거
"""
from fastapi import APIRouter, HTTPException, Depends, UploadFile, File, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from datetime import datetime
from uuid import uuid4
import tempfile
import os

from app.core.config import settings
from app.core.supabase import supabase
from app.services.pdf_parser import pdf_parser_service
from app.services.pubmed_crawler import pubmed_crawler
from app.services.rag_service import rag_service
from app.services.crawler_stealth import AmbeedStealthScraper

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
    payload: Optional[str]
    linker: Optional[str]
    enrichment_source: str
    created_at: str
    outcome_type: Optional[str] = None # Success, Failure, Ongoing
    failure_reason: Optional[str] = None
    is_ai_extracted: bool = False


class ApproveRequest(BaseModel):
    """승인 요청"""
    note: Optional[str] = None


class RejectRequest(BaseModel):
    """반려 요청"""
    reason: str = Field(..., min_length=5, description="반려 사유")


class ManualGoldenSetRequest(BaseModel):
    """골든셋 수동 입력 요청"""
    name: str
    category: str = "antibody" # antibody, payload, linker, clinical_trial
    description: str
    properties: Dict[str, Any]
    outcome_type: Optional[str] = None # Success, Failure, Terminated
    failure_reason: Optional[str] = None
    ip_status: str = "Unknown"


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
    """
    try:
        # 1. Total Users
        users_res = supabase.table("profiles").select("id", count="exact").execute()
        total_users = users_res.count or 0
        
        # 2. Active Simulations (Processing)
        sims_res = supabase.table("projects").select("id", count="exact").eq("status", "processing").execute()
        active_sims = sims_res.count or 0
        
        # 3. Pending Golden Set Items
        pending_res = supabase.table("golden_set_library").select("id", count="exact").eq("status", "draft").execute()
        pending_count = pending_res.count or 0
        
        return {
            "total_users": total_users,
            "active_simulations": active_sims,
            "pending_reviews": pending_count
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/refiner/run")
async def run_ai_refiner(background_tasks: BackgroundTasks):
    """
    AI Refiner 수동 실행 (즉시 트리거)
    """
    try:
        from app.services.ai_refiner import ai_refiner
        # 백그라운드에서 실행 (최대 20개씩 처리)
        background_tasks.add_task(ai_refiner.process_pending_records, max_records=20)
        return {"status": "started", "message": "AI Refiner started in background"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# Credit Management
# ============================================================

@router.post("/credits/grant", response_model=CreditGrantResponse)
async def grant_credits(req: CreditGrantRequest):
    """
    특정 유저 크레딧 지급
    """
    try:
        # 1. Get current balance
        profile_res = supabase.table("profiles").select("credits").eq("id", req.user_id).execute()
        if not profile_res.data:
            raise HTTPException(status_code=404, detail="User not found")
        
        current_credits = profile_res.data[0]["credits"]
        new_balance = current_credits + req.amount
        
        # 2. Update profile
        supabase.table("profiles").update({"credits": new_balance}).eq("id", req.user_id).execute()
        
        # 3. Log transaction
        trans_res = supabase.table("transactions").insert({
            "user_id": req.user_id,
            "amount": req.amount,
            "type": "admin_grant",
            "description": req.reason
        }).execute()
        
        transaction_id = "pending"
        if trans_res.data:
            transaction_id = trans_res.data[0]["id"]

        return {
            "status": "success",
            "user_id": req.user_id,
            "amount": req.amount,
            "new_balance": new_balance,
            "transaction_id": transaction_id
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))# ============================================================
# User Management
# ============================================================

@router.get("/users")
async def get_all_users():
    """
    전체 사용자 목록 조회 (프로필 정보 포함)
    """
    try:
        # profiles 테이블 조회
        try:
            res = supabase.table("profiles").select("*").order("created_at", desc=True).execute()
            return res.data
        except Exception as db_error:
            print(f"Warning: Failed to fetch profiles. Error: {db_error}")
            return []
    except Exception as e:
        print(f"Error in get_all_users: {e}")
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# Golden Set Management (Human-in-the-Loop)
# ============================================================

@router.get("/goldenset/drafts", response_model=List[GoldenSetDraft])
async def get_golden_set_drafts():
    """
    검토 대기 중인 Golden Set 목록 조회
    """
    try:
        response = supabase.table("golden_set_library") \
            .select("*") \
            .eq("status", "draft") \
            .order("created_at", desc=True) \
            .execute()
            
        drafts = []
        for item in response.data:
            drafts.append(GoldenSetDraft(
                id=item.get("id"),
                drug_name=item.get("name") or "Unknown",
                target=item.get("properties", {}).get("target") or "Unknown",
                payload=item.get("properties", {}).get("payload"),
                linker=item.get("properties", {}).get("linker"),
                enrichment_source=item.get("enrichment_source") or "unknown",
                created_at=item.get("created_at"),
                outcome_type=item.get("outcome_type"),
                failure_reason=item.get("failure_reason"),
                is_ai_extracted=item.get("is_ai_extracted") or False
            ))
            
        return drafts
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/goldenset/{golden_set_id}")
async def get_golden_set_detail(golden_set_id: str):
    """
    Golden Set 상세 조회
    """
    try:
        response = supabase.table("golden_set_library") \
            .select("*") \
            .eq("id", golden_set_id) \
            .execute()
            
        if not response.data:
            raise HTTPException(status_code=404, detail="Golden Set not found")
            
        return response.data[0]
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/goldenset/{golden_set_id}/approve")
async def approve_golden_set(golden_set_id: str, req: ApproveRequest, background_tasks: BackgroundTasks):
    """
    Golden Set 승인 및 RAG 인덱싱 (Async)
    """
    try:
        # 1. DB 조회
        item_res = supabase.table("golden_set_library").select("*").eq("id", golden_set_id).execute()
        if not item_res.data:
            raise HTTPException(status_code=404, detail="Golden Set not found")
        
        item = item_res.data[0]
        
        # 2. 상태 업데이트
        update_data = {
            "status": "approved",
            # "approved_at": datetime.utcnow().isoformat() # 컬럼 없음
        }
        
        supabase.table("golden_set_library") \
            .update(update_data) \
            .eq("id", golden_set_id) \
            .execute()
            
        # 3. 비동기 임베딩 생성
        background_tasks.add_task(
            rag_service.index_golden_set_item, 
            golden_set_id, 
            item.get("description", ""), 
            item.get("properties", {})
        )
    
        return {
            "status": "approved",
            "golden_set_id": golden_set_id,
            "message": "Data approved. Indexing started in background."
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/goldenset/{golden_set_id}/reject")
async def reject_golden_set(golden_set_id: str, req: RejectRequest):
    """
    Golden Set 반려
    """
    try:
        update_data = {
            "status": "rejected",
            "reviewer_note": req.reason,
        }
        
        supabase.table("golden_set_library") \
            .update(update_data) \
            .eq("id", golden_set_id) \
            .execute()
            
        return {
            "status": "rejected",
            "golden_set_id": golden_set_id,
            "reason": req.reason
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/goldenset")
async def create_golden_set(req: ManualGoldenSetRequest, background_tasks: BackgroundTasks):
    """
    Golden Set 수동 입력
    """
    try:
        new_entry = {
            "name": req.name,
            "category": req.category,
            "description": req.description,
            "properties": req.properties,
            "status": "approved", # 수동 입력은 바로 승인 처리 (또는 draft로 할 수도 있음)
            "enrichment_source": "manual_entry",
            "outcome_type": req.outcome_type,
            "failure_reason": req.failure_reason,
            "ip_status": req.ip_status
        }
        
        res = supabase.table("golden_set_library").insert(new_entry).execute()
        
        if not res.data:
            raise HTTPException(status_code=500, detail="Failed to insert data")
            
        item_id = res.data[0]['id']
        
        # 바로 인덱싱
        background_tasks.add_task(
            rag_service.index_golden_set_item,
            item_id,
            req.description,
            req.properties
        )
        
        return {
            "status": "success",
            "id": item_id,
            "message": "Golden Set created and indexing started."
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# Prompt Management
# ============================================================

@router.get("/prompts")
async def get_all_prompts():
    """모든 에이전트의 현재 활성 프롬프트 조회"""
    try:
        # is_active=true 인 것만 조회
        res = supabase.table("agent_prompts").select("*").eq("is_active", True).execute()
        return res.data
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/prompts/{agent_id}")
async def get_prompt(agent_id: str):
    """특정 에이전트의 프롬프트 히스토리 조회"""
    try:
        res = supabase.table("agent_prompts") \
            .select("*") \
            .eq("agent_id", agent_id) \
            .order("created_at", desc=True) \
            .execute()
        return res.data
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.put("/prompts/{agent_id}")
async def update_prompt(agent_id: str, req: PromptUpdateRequest):
    """
    에이전트 프롬프트 업데이트 (새 버전 생성)
    """
    try:
        version = req.version or f"v{datetime.now().strftime('%Y%m%d%H%M')}"
        
        # 1. 기존 활성 프롬프트 비활성화 (Optional: 하나만 활성화 유지 정책 시)
        # supabase.table("agent_prompts").update({"is_active": False}).eq("agent_id", agent_id).execute()
        
        # 2. 새 버전 Insert
        new_prompt = {
            "agent_id": agent_id,
            "version": version,
            "system_prompt": req.system_prompt,
            "is_active": True, # 바로 활성화 (또는 draft 로직 분리 가능)
            "is_live": False
        }
        
        res = supabase.table("agent_prompts").insert(new_prompt).execute()
        
        return {
            "status": "success",
            "data": res.data[0]
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# Simulation Logs (User Ops)
# ============================================================

@router.get("/simulations")
async def get_simulation_logs(limit: int = 50):
    """
    전체 시뮬레이션 로그 조회
    """
    try:
        # projects 테이블 조회 + user info join
        # Try joining with profiles first
        try:
            res = supabase.table("projects") \
                .select("*, profiles(email)") \
                .order("created_at", desc=True) \
                .limit(limit) \
                .execute()
            return res.data
        except Exception as join_error:
            print(f"Warning: Failed to join profiles. Fetching projects only. Error: {join_error}")
            # Fallback: Fetch projects without join
            res = supabase.table("projects") \
                .select("*") \
                .order("created_at", desc=True) \
                .limit(limit) \
                .execute()
            return res.data
        
    except Exception as e:
        print(f"Error fetching simulation logs: {e}")
        raise HTTPException(status_code=500, detail=f"Failed to fetch logs: {str(e)}")


# ============================================================
# Health Check
# ============================================================

@router.get("/health")
async def admin_health():
    """관리자 API 헬스 체크"""
    return {"status": "healthy", "service": "admin-api"}


# ============================================================
# Crawler Triggers
# ============================================================

class CrawlerRunRequest(BaseModel):
    category: str = "Payload" # Payload, Linker, Conjugate
    max_pages: int = 1

@router.post("/crawler/ambeed/run")
async def run_ambeed_crawler(req: CrawlerRunRequest, background_tasks: BackgroundTasks):
    """
    Ambeed Stealth Crawler 실행 (Background)
    """
    try:
        scraper = AmbeedStealthScraper()
        background_tasks.add_task(scraper.run, req.category, req.max_pages)
        
        return {
            "status": "started",
            "message": f"Ambeed crawler started for category: {req.category}",
            "mode": "stealth"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/crawler/creative/run")
async def run_creative_crawler(background_tasks: BackgroundTasks, target: str = "HER2"):
    """
    Creative Biolabs Crawler 실행 (Background)
    """
    try:
        from app.services.creative_biolabs_crawler import creative_crawler
        background_tasks.add_task(creative_crawler.run, target)
        
        return {
            "status": "started",
            "message": f"Creative Biolabs crawler started for target: {target}"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# AI Control Dashboard (System Monitoring & Control)
# ============================================================

class SystemStatusRequest(BaseModel):
    """시스템 상태 변경 요청"""
    status: str = Field(..., pattern="^(ACTIVE|PAUSED)$", description="ACTIVE or PAUSED")


@router.get("/system/status")
async def get_system_status():
    """AI Refiner 시스템 상태 조회"""
    try:
        res = supabase.table("system_config").select("value").eq("key", "AI_REFINER_STATUS").execute()
        return {"status": res.data[0]["value"] if res.data else "ACTIVE"}
    except Exception as e:
        # 테이블이 없으면 기본값 반환
        return {"status": "ACTIVE"}


@router.post("/system/status")
async def set_system_status(req: SystemStatusRequest):
    """AI Refiner 시스템 상태 변경 (ACTIVE/PAUSED)"""
    try:
        supabase.table("system_config").upsert({
            "key": "AI_REFINER_STATUS",
            "value": req.status,
            "updated_at": datetime.utcnow().isoformat()
        }).execute()
        return {"status": req.status, "message": f"System status changed to {req.status}"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/refiner/dashboard")
async def get_refiner_dashboard():
    """AI Refiner 대시보드 통합 데이터"""
    try:
        from app.services.cost_tracker import cost_tracker
        from datetime import date
        
        # 1. 시스템 상태
        try:
            status_res = supabase.table("system_config").select("value").eq("key", "AI_REFINER_STATUS").execute()
            system_status = status_res.data[0]["value"] if status_res.data else "ACTIVE"
        except:
            system_status = "ACTIVE"
        
        # 2. 비용 정보
        cost_summary = await cost_tracker.get_usage_summary()
        
        # 3. 큐 상태 - Pending (ai_refined = false)
        pending_res = supabase.table("golden_set_library").select("count", count="exact").eq("ai_refined", False).execute()
        pending_count = pending_res.count or 0
        
        # 4. 오늘 처리된 개수 (ai_refined = true AND today)
        today = date.today().isoformat()
        enriched_res = supabase.table("golden_set_library")\
            .select("count", count="exact")\
            .eq("ai_refined", True)\
            .gte("created_at", f"{today}T00:00:00")\
            .execute()
        enriched_today = enriched_res.count or 0
        
        return {
            "system_status": system_status,
            "cost": cost_summary,
            "queue": {
                "pending": pending_count,
                "enriched_today": enriched_today
            }
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/refiner/recent-logs")
async def get_refiner_recent_logs(limit: int = 5):
    """최근 AI Refiner 처리 결과"""
    try:
        res = supabase.table("golden_set_library")\
            .select("id, name, outcome_type, failure_reason, properties, smiles_code, created_at")\
            .eq("ai_refined", True)\
            .order("created_at", desc=True)\
            .limit(limit)\
            .execute()
        return res.data
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


