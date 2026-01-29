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
import json
import asyncio

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


class AIChatRequest(BaseModel):
    """AI 어시스턴트 채팅 요청"""
    record_id: Any # Relaxed from str to handle int/str
    message: str
    context: Any = None # Relaxed from Dict[str, Any]
    mode: str = "rag" # rag, general
    target_field: Optional[str] = None # For Autofill intent


@router.post("/ai/chat")
async def ai_assistant_chat(req: AIChatRequest):
    """
    AI 어시스턴트 사이드바용 채팅 API (Dual Agent System)
    - RAG Mode: Multi-Source Search + Strict Citations
    - General Mode: LLM Knowledge + Disclaimer
    """
    try:
        import google.generativeai as genai
        from app.services.rag_service import rag_service
        
        if not settings.GOOGLE_API_KEY:
            raise HTTPException(status_code=500, detail="GOOGLE_API_KEY is not configured")
            
        genai.configure(api_key=settings.GOOGLE_API_KEY)
        model_id = settings.GEMINI_MODEL_ID or 'gemini-2.5-flash'
        model = genai.GenerativeModel(model_id)

        # 1. Determine Mode & Context
        mode = req.mode
        context_str = ""
        
        # Autofill Override
        if req.target_field:
            mode = "rag" # Autofill always uses RAG
            req.message = f"Find the value for the field '{req.target_field}' of this item. If found, provide ONLY the value and the source."

        if mode == "rag":
            # Multi-Source Search
            search_query = req.message
            # If record context is provided, enrich query with drug name
            if req.context:
                drug_name = req.context.get("drug_name") or req.context.get("product_name") or req.context.get("name")
                if drug_name:
                    search_query = f"{drug_name} {req.message}"

            retrieved_chunks = await rag_service.search_all(search_query, top_k_per_source=3)

            # Format Context with Source IDs
            context_str = "RETRIEVED KNOWLEDGE:\n"
            for chunk in retrieved_chunks:
                source_id = chunk.get('id')
                source_type = chunk.get('source', 'unknown')
                content = chunk.get('content', '')
                metadata = chunk.get('metadata', {})
                context_str += f"- [{source_type}:{source_id}] {content} (Meta: {metadata})\n"

            system_prompt = """You are a STRICT Data Analyst and Scientific Literature Expert.
Answer the user's question using the provided RETRIEVED KNOWLEDGE and your scientific expertise.

CRITICAL CITATION RULES (MUST FOLLOW):
1. Every factual statement MUST include a citation.
2. For internal data: Use format `[Source: table:id]`
3. For scientific literature: ALWAYS include PMID when available.
   - Format: `[PMID: 12345678]` or `[DOI: 10.xxxx/xxxxx]`
   - If you mention a study or paper, you MUST provide the PMID link: https://pubmed.ncbi.nlm.nih.gov/PMID/
4. If no PMID is available but you reference general knowledge, explicitly state: "[No PMID - General Knowledge]"

SMILES VALIDATION RULES:
- When providing SMILES structures, always format as: `SMILES: [structure]`
- Include molecular weight (MW) when discussing chemical structures
- Mention if the SMILES is from PubChem, ChEMBL, or other databases

EXAMPLE RESPONSE FORMAT:
"Trastuzumab (Herceptin) is a humanized monoclonal antibody targeting HER2/neu receptor [PMID: 9851916](https://pubmed.ncbi.nlm.nih.gov/9851916/).
The recommended SMILES for MMAE is: SMILES: CC(C)C[C@H]... [Source: commercial_reagents:abc123]"

If the answer is not in the context, say "내부 데이터에서 찾을 수 없습니다. 일반 지식 모드로 전환해주세요."
"""
        else:
            # General Mode - Still require PMID citations
            system_prompt = """You are a Helpful AI Assistant specialized in ADCs and pharmaceuticals.
Answer using your general knowledge but ALWAYS provide academic citations.

**[General Knowledge Mode]** - This answer is based on general AI knowledge, not internal data.

CRITICAL CITATION RULES (MUST FOLLOW):
1. For ANY scientific claim, provide a PMID or DOI reference.
2. Format: [PMID: 12345678](https://pubmed.ncbi.nlm.nih.gov/12345678/)
3. If you cannot provide a specific PMID, indicate: "[Citation needed - verify independently]"
4. When discussing chemical structures, provide SMILES in format: `SMILES: [structure]`

Your response should be:
- Scientifically accurate
- Include at least 1-3 relevant PMID citations
- Clearly formatted with structure information when applicable
"""

        # 2. Construct Full Prompt
        full_prompt = f"""{system_prompt}

USER CONTEXT (Current Item):
{json.dumps(req.context or {}, indent=2, ensure_ascii=False)[:2000]}

{context_str}

User Question: {req.message}
"""
        
        # 3. Call Gemini
        try:
            response = await asyncio.get_event_loop().run_in_executor(
                None, lambda: model.generate_content(full_prompt)
            )
            answer = response.text.strip() if response and response.text else "No response generated."
            
            return {
                "answer": answer,
                "mode": mode,
                "sources": retrieved_chunks if mode == "rag" else []
            }
            
        except Exception as e:
            print(f"Gemini Error: {e}")
            return {"answer": "Error generating response."}

    except Exception as e:
        print(f"AI Chat Critical Error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/stats")

@router.get("/stats")
async def get_admin_stats():
    """
    대시보드 KPI 조회 (Expanded)
    """
    try:
        # 1. Design Engine Stats
        design_res = supabase.table("design_runs").select("status", count="exact").execute()
        design_stats = {
            "pending": 0,
            "processing": 0,
            "completed": 0
        }
        if design_res.data:
            for item in design_res.data:
                status = item.get("status", "pending")
                design_stats[status] = design_stats.get(status, 0) + 1

        # 2. Data Maturity Stats
        # Golden Set
        gs_approved = supabase.table("golden_set_library").select("id", count="exact").eq("status", "approved").execute()
        gs_draft = supabase.table("golden_set_library").select("id", count="exact").eq("status", "draft").execute()
        
        # Knowledge Base
        kb_total = supabase.table("knowledge_base").select("id", count="exact").execute()
        kb_completed = supabase.table("knowledge_base").select("id", count="exact").eq("rag_status", "completed").execute()

        data_stats = {
            "golden_set_approved": gs_approved.count or 0,
            "golden_set_candidates": gs_draft.count or 0,
            "kb_total": kb_total.count or 0,
            "kb_completed": kb_completed.count or 0
        }

        # 3. Inventory Stats
        antibodies = supabase.table("antibody_library").select("id", count="exact").execute()
        reagents = supabase.table("commercial_reagents").select("id", count="exact").execute()
        targets = supabase.table("target_master").select("id", count="exact").execute()

        inventory_stats = {
            "antibodies": antibodies.count or 0,
            "reagents": reagents.count or 0,
            "targets": targets.count or 0
        }

        # 4. System Stats
        # AI Cost (Mock calculation based on logs if available, or just sum usage)
        # For now, we'll try to sum 'cost' column if it exists, or just count logs
        # Assuming llm_usage_logs has a 'cost' column. If not, we mock it.
        ai_cost_month = 0.0
        try:
            # This is a placeholder. In real app, sum(cost) where created_at > 1 month ago
            # Supabase JS client doesn't support sum easily without RPC. 
            # We will just count tokens or use a mock value for now if RPC not ready.
            logs_res = supabase.table("llm_usage_logs").select("total_tokens").order("created_at", desc=True).limit(100).execute()
            if logs_res.data:
                # Rough estimate: $0.00001 per token
                total_tokens = sum(item.get('total_tokens', 0) for item in logs_res.data)
                ai_cost_month = total_tokens * 0.00001
        except:
            ai_cost_month = 45.20 # Fallback mock

        # Data Velocity (24h)
        one_day_ago = (datetime.utcnow().replace(microsecond=0) ).isoformat() # Simplified
        # Actually need proper timedelta, importing timedelta
        from datetime import timedelta
        one_day_ago = (datetime.utcnow() - timedelta(days=1)).isoformat()
        
        velocity_res = supabase.table("knowledge_base").select("id", count="exact").gte("created_at", one_day_ago).execute()
        data_velocity_24h = velocity_res.count or 0

        system_stats = {
            "ai_cost_month": round(ai_cost_month, 2),
            "data_velocity_24h": data_velocity_24h
        }

        # 5. Recent Activity (Live Feed)
        # Fetch recent KB items
        kb_recent = supabase.table("knowledge_base").select("id, title, created_at, source_type").order("created_at", desc=True).limit(3).execute()
        # Fetch recent Golden Set items
        gs_recent = supabase.table("golden_set_library").select("id, name, created_at, status").order("created_at", desc=True).limit(3).execute()

        recent_activity = []
        for item in (kb_recent.data or []):
            recent_activity.append({
                "type": "paper",
                "id": item.get("id"),
                "title": item.get("title"),
                "time": item.get("created_at"),
                "desc": f"New {item.get('source_type')} article collected"
            })
        
        for item in (gs_recent.data or []):
            recent_activity.append({
                "type": "goldenset",
                "id": item.get("id"),
                "title": item.get("name"),
                "time": item.get("created_at"),
                "desc": f"Golden Set {item.get('status')}"
            })
        
        # Sort by time desc and take top 5
        recent_activity.sort(key=lambda x: x['time'], reverse=True)
        recent_activity = recent_activity[:5]

        return {
            "design_stats": design_stats,
            "data_stats": data_stats,
            "inventory_stats": inventory_stats,
            "system_stats": system_stats,
            "recent_activity": recent_activity
        }
    except Exception as e:
        print(f"Stats Error: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/refiner/run")
async def run_ai_refiner(
    background_tasks: BackgroundTasks, 
    limit: int = 50, 
    mode: str = "partial", # partial, full, daily_import
    source: Optional[str] = None # clinical_trials, open_fda_api
):
    """
    AI Refiner 수동 실행 (즉시 트리거)
    - mode="partial": limit 만큼 실행 (기본 50)
    - mode="full": 대기 중인 모든 항목 실행 (Safety Net: Max 10,000)
    - mode="daily_import": Daily Bulk Import 실행
    """
    try:
        from app.services.ai_refiner import ai_refiner
        
        # Commercial Reagents Batch Logic
        if source == "commercial":
            from app.services.commercial_refiner_service import batch_refine_commercial_reagents
            
            # Count pending
            query = supabase.table("commercial_reagents").select("count", count="exact").eq("ai_refined", False)
            pending_res = query.execute()
            pending_count = pending_res.count or 0
            
            target_limit = limit
            if mode == "full":
                target_limit = min(pending_count, 10000)
            
            # Run in background
            background_tasks.add_task(batch_refine_commercial_reagents, limit=target_limit, mode=mode)
            
            return {
                "status": "started", 
                "message": f"Commercial Reagents Refiner started ({target_limit} items)",
                "count": target_limit,
                "mode": mode
            }

        if mode == "daily_import":
            # Source에 따라 Import 로직 분기
            if source == "open_fda_api":
                from app.services.openfda_service import openfda_service
                job_id = f"sync_openfda_{uuid4().hex[:8]}"
                # DB에 잡 등록
                supabase.table("sync_jobs").insert({
                    "id": job_id, "status": "queued", "source": "openfda", "started_at": datetime.utcnow().isoformat()
                }).execute()
                background_tasks.add_task(openfda_service.sync_to_db, job_id, mode="daily", limit=100)
                return {
                    "status": "started", 
                    "message": "OpenFDA Daily Import started in background", 
                    "mode": "daily_import",
                    "source": "open_fda_api"
                }
            else:
                from app.services.bulk_importer import BulkImporter
                importer = BulkImporter()
                background_tasks.add_task(importer.run_import, max_studies=5000, mode="daily")
                return {
                    "status": "started", 
                    "message": "ClinicalTrials Daily Import started in background",
                    "mode": "daily_import"
                }
            
        # Full 모드일 경우 pending 개수 조회
        target_limit = limit
        if mode == "full":
            # Pending 개수 조회
            query = supabase.table("golden_set_library").select("count", count="exact").eq("ai_refined", False)
            if source:
                query = query.eq("enrichment_source", source)
            pending_res = query.execute()
            pending_count = pending_res.count or 0
            
            if pending_count == 0:
                return {"status": "skipped", "message": "No pending items to refine"}
                
            # Safety Net: 최대 10,000개
            target_limit = min(pending_count, 10000)
            
        # 백그라운드에서 실행
        background_tasks.add_task(ai_refiner.process_pending_records, max_records=target_limit, source_filter=source)
        
        # 예상 소요 시간 (개당 2초)
        estimated_seconds = target_limit * 2
        
        return {
            "status": "started", 
            "message": f"AI Refiner started in background ({target_limit} items, source: {source or 'All'})",
            "count": target_limit,
            "estimated_seconds": estimated_seconds,
            "mode": mode
        }
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

@router.get("/goldenset/library")
async def get_golden_set_library(
    limit: int = 20,
    offset: int = 0,
    tab: str = "success", # success, failure
    search: str = ""
):
    """
    골든셋 라이브러리 (사장님이 승격시킨 데이터) 조회
    - tab='success': status='approved' AND outcome_type='Success'
    - tab='failure': status='approved' AND (outcome_type!='Success' OR outcome_type IS NULL)
    """
    try:
        query = supabase.table("golden_set_library").select("*", count="exact")
        
        # 기본 조건: 승격된 데이터 (status='approved')
        # (만약 '승격'이라는 별도 플래그가 없다면 status로 구분)
        query = query.eq("status", "approved")
        
        if tab == "success":
            query = query.eq("outcome_type", "Success")
        else: # failure tab
            # outcome_type이 Success가 아닌 모든 것 (Failure, Unknown, Ongoing 등)
            query = query.neq("outcome_type", "Success")
            
        if search:
            query = query.ilike("name", f"%{search}%")
            
        response = query.order("updated_at", desc=True).range(offset, offset + limit - 1).execute()
        
        return {
            "data": response.data,
            "total": response.count,
            "limit": limit,
            "offset": offset
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/goldenset/drafts")
async def get_golden_set_drafts(
    limit: int = 20,
    offset: int = 0,
    search: str = "",
    source: str = ""
):
    """
    검토 대기 중인 Golden Set 목록 조회 (페이지네이션 + 검색/필터)
    - limit: 한 페이지당 개수 (기본 20)
    - offset: 오프셋 (페이지 * limit)
    - search: 약물명 검색 (ILIKE)
    - source: enrichment_source 필터 (open_fda_api, clinical_trials_api_v2 등)
    """
    try:
        # Count 쿼리 (전체 개수)
        count_query = supabase.table("golden_set_library") \
            .select("*", count="exact") \
            .eq("status", "draft")
        
        if search:
            count_query = count_query.ilike("name", f"%{search}%")
        if source:
            count_query = count_query.eq("enrichment_source", source)
        
        count_result = count_query.execute()
        total_count = count_result.count or 0
        
        # 데이터 쿼리 (페이지네이션)
        data_query = supabase.table("golden_set_library") \
            .select("*") \
            .eq("status", "draft")
        
        if search:
            data_query = data_query.ilike("name", f"%{search}%")
        if source:
            data_query = data_query.eq("enrichment_source", source)
        
        response = data_query \
            .order("created_at", desc=True) \
            .range(offset, offset + limit - 1) \
            .execute()
            
        drafts = []
        for item in response.data:
            try:
                drafts.append({
                    "id": item.get("id"),
                    "drug_name": item.get("name") or "Unknown",
                    "target": item.get("properties", {}).get("target") or "Unknown",
                    "payload": item.get("properties", {}).get("payload"),
                    "linker": item.get("properties", {}).get("linker"),
                    "smiles_code": item.get("smiles_code"),  # SMILES 추가
                    "enrichment_source": item.get("enrichment_source") or "unknown",
                    "created_at": item.get("created_at"),
                    "outcome_type": item.get("outcome_type"),
                    "failure_reason": item.get("failure_reason"),
                    "binding_affinity": item.get("binding_affinity") or item.get("properties", {}).get("binding_affinity"),
                    "isotype": item.get("isotype") or item.get("properties", {}).get("isotype"),
                    "host_species": item.get("host_species") or item.get("properties", {}).get("host_species"),
                    "orr_pct": item.get("orr_pct") or item.get("properties", {}).get("orr_pct"),
                    "os_months": item.get("os_months") or item.get("properties", {}).get("os_months"),
                    "pfs_months": item.get("pfs_months") or item.get("properties", {}).get("pfs_months"),
                    "is_ai_extracted": item.get("ai_refined") or False,
                    "raw_data": item.get("raw_data", {}),
                    "properties": item.get("properties", {}) or {}  # 전체 properties도 포함, None일 경우 빈 dict
                })
            except Exception as item_error:
                print(f"Skipping bad item {item.get('id', 'unknown')}: {item_error}")
                continue
            
        return {
            "data": drafts,
            "total": total_count,
            "limit": limit,
            "offset": offset
        }
        
    except Exception as e:
        print(f"Critical Error in get_golden_set_drafts: {e}")
        raise HTTPException(status_code=500, detail=f"Server Error: {str(e)}")


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


@router.post("/goldenset/{golden_set_id}/refine")
async def refine_single_golden_set(golden_set_id: str, background_tasks: BackgroundTasks):
    """
    개별 레코드 AI 재분석 + SMILES 생성
    스테이징 에어리어에서 즉시 분석 버튼용
    """
    try:
        from app.services.ai_refiner import ai_refiner
        
        # 1. 레코드 조회
        item_res = supabase.table("golden_set_library").select("*").eq("id", golden_set_id).execute()
        if not item_res.data:
            raise HTTPException(status_code=404, detail="Record not found")
        
        record = item_res.data[0]
        
        # 2. AI Refiner로 분석 (동기 실행)
        analysis = await ai_refiner.refine_single_record(record)
        
        if analysis and "error" not in analysis:
            drug_name = analysis.get("drug_name") or record.get("name")
            
            # 3. SMILES 생성 (PubChem → AI Fallback)
            smiles_data = {}
            if drug_name and not record.get("smiles_code"):
                smiles_data = await ai_refiner.enrich_with_pubchem(drug_name)
                if not smiles_data or "error" in smiles_data:
                    smiles_data = await ai_refiner.generate_smiles_with_ai(drug_name)
            
            # 4. DB 업데이트
            existing_props = record.get("properties", {}) or {}
            existing_props["ai_analysis"] = analysis
            
            update_payload = {
                "name": drug_name or record.get("name"),
                "outcome_type": analysis.get("outcome_type", "Unknown"),
                "failure_reason": analysis.get("failure_reason"),
                "relevance_score": analysis.get("relevance_score", 0.0),
                "ai_refined": True,
                "properties": existing_props
            }
            
            # SMILES 추가
            if smiles_data and "error" not in smiles_data:
                update_payload["smiles_code"] = smiles_data.get("smiles_code")
            
            supabase.table("golden_set_library").update(update_payload).eq("id", golden_set_id).execute()
            
            return {
                "status": "success",
                "id": golden_set_id,
                "analysis": analysis,
                "smiles_code": update_payload.get("smiles_code") or record.get("smiles_code"),
                "message": "AI 분석 완료"
            }
        else:
            error_msg = analysis.get("error") if analysis else "Unknown Error"
            return {
                "status": "error",
                "id": golden_set_id,
                "message": f"AI 분석 실패: {error_msg}"
            }
        
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
    Uses the advanced AmbeedCrawler with AI enrichment.
    """
    try:
        from app.services.ambeed_crawler import ambeed_crawler
        
        # Map frontend category to crawler search term if needed, or pass directly
        # Frontend sends: "Payload", "Linker", "Conjugate"
        # Crawler expects keys in CATEGORIES or "all"
        
        # Mapping logic (can be moved to crawler, but doing here for clarity)
        category_map = {
            "Payload": "ADC Toxins", # or ADC Cytotoxin
            "Linker": "ADC Linkers",
            "Conjugate": "ADC Cytotoxin" # Placeholder if no specific conjugate category
        }
        
        search_term = category_map.get(req.category, req.category)
        
        # Convert max_pages to limit (approx 20 items per page)
        limit = req.max_pages * 20
        
        job_id = f"crawl_ambeed_{uuid4().hex[:8]}"
        
        data = {
            "id": job_id, 
            "status": "queued", 
            "source": "ambeed", 
            "started_at": datetime.utcnow().isoformat()
        }
        supabase.table("sync_jobs").insert(data).execute()
        
        background_tasks.add_task(ambeed_crawler.run, search_term, limit, job_id)
        
        return {
            "status": "started",
            "message": f"Ambeed crawler started for {req.category} (limit={limit})",
            "job_id": job_id,
            "mode": "stealth_ai"
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
        
        # 3.1 OpenFDA Pending
        openfda_pending_res = supabase.table("golden_set_library").select("count", count="exact").eq("ai_refined", False).eq("enrichment_source", "open_fda_api").execute()
        openfda_pending = openfda_pending_res.count or 0
        
        # 3.2 Commercial Pending
        comm_pending_res = supabase.table("commercial_reagents").select("count", count="exact").eq("ai_refined", False).execute()
        comm_pending = comm_pending_res.count or 0
        
        # 4. 오늘 처리된 개수 (ai_refined = true AND today)
        today = date.today().isoformat()
        enriched_res = supabase.table("golden_set_library")\
            .select("count", count="exact")\
            .eq("ai_refined", True)\
            .gte("updated_at", f"{today}T00:00:00")\
            .execute()
        enriched_today = enriched_res.count or 0
        
        return {
            "system_status": system_status,
            "cost": cost_summary,
            "queue": {
                "pending": pending_count,
                "openfda_pending": openfda_pending,
                "commercial_pending": comm_pending,
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


# ============================================================
# Knowledge Base Refiner (PubMed/BioRxiv)
# ============================================================

@router.post("/knowledge/refine")
async def run_knowledge_refiner(background_tasks: BackgroundTasks, limit: int = 20):
    """
    Knowledge Base Refiner 수동 실행 (PubMed/BioRxiv 분석)
    """
    try:
        from app.services.knowledge_refiner import knowledge_refiner
        background_tasks.add_task(knowledge_refiner.process_pending_items, batch_size=limit)
        return {"status": "started", "message": f"Knowledge Refiner started (limit: {limit})"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/knowledge/logs")
async def get_knowledge_logs(limit: int = 10):
    """최근 Knowledge Base Refiner 처리 결과 (Spot Check용)"""
    try:
        # rag_status가 processed인 항목 조회
        res = supabase.table("knowledge_base")\
            .select("id, title, summary, relevance_score, ai_reasoning, updated_at")\
            .eq("rag_status", "processed")\
            .order("updated_at", desc=True)\
            .limit(limit)\
            .execute()
            
        # 프론트엔드 포맷에 맞게 매핑 (RecentLog 인터페이스 호환)
        logs = []
        for item in res.data:
            logs.append({
                "id": item.get("id"),
                "name": item.get("title")[:50] + "..." if item.get("title") else "Unknown",
                "outcome_type": f"Score: {item.get('relevance_score', 0)}",
                "failure_reason": item.get("ai_reasoning"), # Spot Check 모달에서 Reason으로 표시됨
                "properties": {},
                "smiles_code": None,
                "created_at": item.get("updated_at")
            })
            
        return logs
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# Total Data Inventory Management (New Feature)
# ============================================================

@router.get("/inventory/commercial")
async def get_commercial_inventory(
    limit: int = 20,
    offset: int = 0,
    search: str = "",
    missing_data_only: bool = False
):
    """
    상용 시약 전체 조회 (필터링 포함)
    """
    try:
        query = supabase.table("commercial_reagents").select("*", count="exact")
        
        if search:
            query = query.ilike("product_name", f"%{search}%")
            
        if missing_data_only:
            # 타겟이 없거나 SMILES가 없는 경우
            # Supabase-py OR filter logic is tricky, usually raw string:
            # .or_("target.is.null,smiles_code.is.null")
            query = query.or_("target.is.null,smiles_code.is.null")
            
        response = query.order("crawled_at", desc=True).range(offset, offset + limit - 1).execute()
        
        return {
            "data": response.data,
            "total": response.count,
            "limit": limit,
            "offset": offset
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/inventory/knowledge")
async def get_knowledge_inventory(
    limit: int = 20,
    offset: int = 0,
    search: str = "",
    missing_data_only: bool = False
):
    """
    지식 베이스 전체 조회
    """
    try:
        query = supabase.table("knowledge_base").select("*", count="exact")
        
        if search:
            query = query.ilike("title", f"%{search}%")
            
        if missing_data_only:
             # 요약이나 추론 결과가 없는 경우
            query = query.or_("summary.is.null,ai_reasoning.is.null")
            
        response = query.order("created_at", desc=True).range(offset, offset + limit - 1).execute()
        
        return {
            "data": response.data,
            "total": response.count,
            "limit": limit,
            "offset": offset
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

class PromoteRequest(BaseModel):
    target_table: str = "golden_set_library" # Default
    note: Optional[str] = None

@router.post("/inventory/promote/{id}")
async def promote_to_golden_set(id: str, req: PromoteRequest):
    """
    상용 시약 데이터를 골든셋으로 승격 (Copy & Insert)
    """
    try:
        # 1. Fetch source data
        src_res = supabase.table("commercial_reagents").select("*").eq("id", id).execute()
        if not src_res.data:
            raise HTTPException(status_code=404, detail="Item not found")
        
        item = src_res.data[0]
        
        # 2. Map to Golden Set Schema
        new_entry = {
            "name": item.get("product_name"),
            "category": "reagent", # or map from item['category']
            "description": item.get("summary") or item.get("body_text") or "Promoted from commercial reagents",
            "properties": item.get("properties", {}),
            "status": "approved", # Auto-approve promoted items? Or "draft"? User implied "Direct promote"
            "enrichment_source": f"promoted_from_{item.get('source_name', 'commercial')}",
            "outcome_type": "Success", # Default assumption for promoted reagents? Or Unknown
            "target": item.get("target"),
            "smiles_code": item.get("smiles_code"),
            "ai_refined": True
        }
        
        # 3. Insert into Golden Set
        res = supabase.table("golden_set_library").insert(new_entry).execute()
        
        # 4. Optional: Mark original as promoted? (Not strictly required by schema)
        
        return {"status": "success", "new_id": res.data[0]['id'], "message": "Item promoted to Golden Set"}
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

class PatchRequest(BaseModel):
    updates: Dict[str, Any]

@router.patch("/inventory/{table}/{id}")
async def patch_inventory_item(table: str, id: str, req: PatchRequest, background_tasks: BackgroundTasks):
    """
    공통 인라인 수정 (수동 편집)
    수정 시 자동으로 Embedding 재생성 및 수동 수정 플래그 설정
    """
    valid_tables = ["commercial_reagents", "knowledge_base", "golden_set_library"]
    if table not in valid_tables:
        raise HTTPException(status_code=400, detail="Invalid table")
        
    try:
        updates = req.updates
        
        # 0. Field Mapping (Handle missing columns by putting into properties)
        # These fields are not columns in the table but should be stored in the 'properties' JSONB column
        property_fields = [
            # Bio Metrics
            "binding_affinity", "isotype", "host_species", "orr_pct", "os_months", "pfs_months", 
            "dor_months", "patient_count", "adverse_events_grade3_pct",
            # ADC Design
            "target_1", "target_2", "target_symbol", "antibody_format", "linker_type", "dar", "gene_id", "uniprot_id",
            # Chemical / IP
            "payload_smiles", "linker_smiles", "full_smiles", "canonical_smiles", "molecular_weight", "patent_id", "patent_expiry"
        ]
        
        # 1. 수동 수정 플래그 설정
        if table in ["commercial_reagents", "golden_set_library"]:
            updates["is_manual_override"] = True
            updates["ai_refined"] = True  # 사람이 고쳤으므로 정제된 것으로 간주
            
            # golden_set_library에서 아직 컬럼이 없는 경우 properties로 이동
            if table == "golden_set_library":
                # 기존 레코드의 properties 가져오기
                existing = supabase.table(table).select("properties").eq("id", id).execute()
                props = (existing.data[0].get("properties") or {}) if existing.data else {}
                
                changed = False
                for field in property_fields:
                    if field in updates:
                        props[field] = updates[field]
                        # columns don't exist yet, so we remove from top-level updates to avoid PGRST204
                        del updates[field]
                        changed = True
                
                # Handle nested properties if sent as a dict (just in case)
                if "properties" in updates and isinstance(updates["properties"], dict):
                    props.update(updates["properties"])
                    del updates["properties"]
                    changed = True

                if changed:
                    updates["properties"] = props

        # 2. 임베딩 실시간 재생성 (RAG 동기화)
        # 수정된 내용 + 기존 내용을 합쳐서 새로운 벡터 생성
        try:
            # 기존 레코드 조회 (텍스트 조합용)
            existing = supabase.table(table).select("*").eq("id", id).execute()
            if existing.data:
                item = existing.data[0]
                # 새로운 텍스트 조합 (수정된 값 우선, 없으면 기존 값)
                name = updates.get("product_name") or updates.get("name") or updates.get("title") or item.get("product_name") or item.get("name") or item.get("title") or ""
                target = updates.get("target") or item.get("target") or ""
                summary = updates.get("summary") or item.get("summary") or ""
                
                # 정량 지표 추가 (검색 품질 향상)
                kd = updates.get("binding_affinity") or item.get("binding_affinity") or ""
                
                combined_text = f"{name} {target} {summary} {kd}".strip()
                
                if combined_text:
                    new_embedding = await rag_service.generate_embedding(combined_text)
                    if new_embedding:
                        updates["embedding"] = new_embedding
        except Exception as embed_error:
            print(f"Embedding regeneration failed: {embed_error}")
            # 임베딩 실패해도 데이터 저장은 진행

        # 3. DB 업데이트
        res = supabase.table(table).update(updates).eq("id", id).execute()
        
        return {"status": "success", "message": "Item updated and re-indexed", "data": res.data[0] if res.data else None}

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# v2.2 Management - Engine Control & Emergency Brake
# ============================================================

class EngineConfigUpdate(BaseModel):
    """엔진 설정 업데이트"""
    useNimApi: Optional[bool] = None
    fallbackEnabled: Optional[bool] = None


class EmergencyStopRequest(BaseModel):
    """긴급 정지 요청"""
    activate: bool


class BudgetConfigUpdate(BaseModel):
    """예산 설정 업데이트"""
    provider: str
    dailyLimitUsd: Optional[float] = None
    monthlyLimitUsd: Optional[float] = None
    alertThresholdPct: Optional[int] = None
    autoFallbackEnabled: Optional[bool] = None
    alertEmail: Optional[str] = None


@router.get("/engine/config")
async def get_engine_config():
    """AI 엔진 설정 조회"""
    try:
        result = {
            "useNimApi": True,
            "fallbackEnabled": True,
            "emergencyStop": False
        }

        # system_config 테이블에서 조회
        configs = ["USE_NIM_API", "FALLBACK_ENABLED", "EMERGENCY_STOP"]
        for key in configs:
            try:
                res = supabase.table("system_config").select("value").eq("key", key).execute()
                if res.data:
                    val = res.data[0]["value"]
                    if key == "USE_NIM_API":
                        result["useNimApi"] = val == "true" or val == True
                    elif key == "FALLBACK_ENABLED":
                        result["fallbackEnabled"] = val == "true" or val == True
                    elif key == "EMERGENCY_STOP":
                        result["emergencyStop"] = val == "true" or val == True
            except:
                continue

        return result
    except Exception as e:
        print(f"Engine config error: {e}")
        return {
            "useNimApi": True,
            "fallbackEnabled": True,
            "emergencyStop": False
        }


@router.post("/engine/switch")
async def switch_engine(req: EngineConfigUpdate):
    """AI 엔진 전환"""
    try:
        if req.useNimApi is not None:
            supabase.table("system_config").upsert({
                "key": "USE_NIM_API",
                "value": "true" if req.useNimApi else "false",
                "updated_at": datetime.utcnow().isoformat()
            }).execute()

            # Audit log
            try:
                supabase.table("admin_audit_logs").insert({
                    "action_type": "ENGINE_SWITCH",
                    "action_category": "system",
                    "target_type": "system_config",
                    "target_id": "USE_NIM_API",
                    "after_value": {"useNimApi": req.useNimApi}
                }).execute()
            except:
                pass

        return {"status": "success", "useNimApi": req.useNimApi}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/budget/emergency-stop")
async def emergency_stop(req: EmergencyStopRequest):
    """긴급 정지 활성화/비활성화"""
    try:
        supabase.table("system_config").upsert({
            "key": "EMERGENCY_STOP",
            "value": "true" if req.activate else "false",
            "updated_at": datetime.utcnow().isoformat()
        }).execute()

        # Audit log
        try:
            supabase.table("admin_audit_logs").insert({
                "action_type": "EMERGENCY_STOP",
                "action_category": "system",
                "target_type": "system_config",
                "target_id": "EMERGENCY_STOP",
                "after_value": {"activated": req.activate}
            }).execute()
        except:
            pass

        return {
            "status": "success",
            "emergencyStop": req.activate,
            "message": "Emergency stop activated" if req.activate else "Emergency stop deactivated"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/budget/status")
async def get_budget_status():
    """예산 상태 조회"""
    try:
        res = supabase.table("budget_configs").select("*").execute()

        budgets = []
        for item in (res.data or []):
            budgets.append({
                "apiProvider": item.get("api_provider"),
                "dailyLimitUsd": float(item.get("daily_limit_usd", 0)),
                "currentUsageUsd": float(item.get("current_usage_usd", 0)),
                "monthlyLimitUsd": float(item.get("monthly_limit_usd", 0)),
                "currentMonthlyUsd": float(item.get("current_monthly_usd", 0)),
                "autoFallbackEnabled": item.get("auto_fallback_enabled", True),
                "alertThresholdPct": item.get("alert_threshold_pct", 80),
                "alertEmail": item.get("alert_email", "")
            })

        # Default budgets if none exist
        if not budgets:
            budgets = [
                {
                    "apiProvider": "nvidia_nim",
                    "dailyLimitUsd": 500,
                    "currentUsageUsd": 0,
                    "monthlyLimitUsd": 10000,
                    "currentMonthlyUsd": 0,
                    "autoFallbackEnabled": True,
                    "alertThresholdPct": 80,
                    "alertEmail": ""
                },
                {
                    "apiProvider": "gemini",
                    "dailyLimitUsd": 200,
                    "currentUsageUsd": 0,
                    "monthlyLimitUsd": 4000,
                    "currentMonthlyUsd": 0,
                    "autoFallbackEnabled": True,
                    "alertThresholdPct": 80,
                    "alertEmail": ""
                }
            ]

        return {"budgets": budgets}
    except Exception as e:
        print(f"Budget status error: {e}")
        return {"budgets": []}


@router.put("/budget/config")
async def update_budget_config(req: BudgetConfigUpdate):
    """예산 설정 업데이트"""
    try:
        updates = {}
        if req.dailyLimitUsd is not None:
            updates["daily_limit_usd"] = req.dailyLimitUsd
        if req.monthlyLimitUsd is not None:
            updates["monthly_limit_usd"] = req.monthlyLimitUsd
        if req.alertThresholdPct is not None:
            updates["alert_threshold_pct"] = req.alertThresholdPct
        if req.autoFallbackEnabled is not None:
            updates["auto_fallback_enabled"] = req.autoFallbackEnabled
        if req.alertEmail is not None:
            updates["alert_email"] = req.alertEmail

        updates["updated_at"] = datetime.utcnow().isoformat()

        supabase.table("budget_configs").update(updates).eq("api_provider", req.provider).execute()

        # Audit log
        try:
            supabase.table("admin_audit_logs").insert({
                "action_type": "BUDGET_UPDATE",
                "action_category": "system",
                "target_type": "budget_configs",
                "target_id": req.provider,
                "after_value": updates
            }).execute()
        except:
            pass

        return {"status": "success", "message": f"Budget config updated for {req.provider}"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


# ============================================================
# v2.2 Management - Refinement Hub (Quarantine & Batch AI Fixer)
# ============================================================

@router.get("/refinement/quarantined")
async def get_quarantined_data(
    status: str = "all",
    source: str = "all",
    limit: int = 50,
    offset: int = 0
):
    """격리 데이터 목록 조회"""
    try:
        query = supabase.table("quarantined_data").select("*", count="exact")

        if status != "all":
            query = query.eq("status", status)
        if source != "all":
            query = query.eq("source_table", source)

        res = query.order("created_at", desc=True).range(offset, offset + limit - 1).execute()

        items = []
        for item in (res.data or []):
            items.append({
                "id": item.get("id"),
                "sourceTable": item.get("source_table"),
                "sourceId": item.get("source_id"),
                "errorType": item.get("error_type"),
                "errorDetails": item.get("error_details"),
                "originalData": item.get("original_data", {}),
                "suggestedFix": item.get("suggested_fix"),
                "status": item.get("status", "pending"),
                "priority": item.get("priority", "normal"),
                "createdAt": item.get("created_at")
            })

        return {"items": items, "total": res.count or 0}
    except Exception as e:
        print(f"Quarantine fetch error: {e}")
        return {"items": [], "total": 0}


@router.post("/refinement/quarantined/{item_id}/approve")
async def approve_quarantined_item(item_id: str):
    """격리 항목 수정 승인 및 적용"""
    try:
        # 1. Get quarantined item
        res = supabase.table("quarantined_data").select("*").eq("id", item_id).execute()
        if not res.data:
            raise HTTPException(status_code=404, detail="Item not found")

        item = res.data[0]
        source_table = item.get("source_table")
        source_id = item.get("source_id")
        suggested_fix = item.get("suggested_fix")

        # 2. Apply fix to source table
        if suggested_fix and source_id:
            supabase.table(source_table).update(suggested_fix).eq("id", source_id).execute()

        # 3. Update quarantined status
        supabase.table("quarantined_data").update({
            "status": "approved",
            "reviewed_at": datetime.utcnow().isoformat()
        }).eq("id", item_id).execute()

        # Audit log
        try:
            supabase.table("admin_audit_logs").insert({
                "action_type": "DATA_APPROVE",
                "action_category": "data",
                "target_type": "quarantined_data",
                "target_id": item_id,
                "before_value": item.get("original_data"),
                "after_value": suggested_fix
            }).execute()
        except:
            pass

        return {"status": "success", "message": "Fix approved and applied"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/refinement/quarantined/{item_id}/reject")
async def reject_quarantined_item(item_id: str):
    """격리 항목 거부"""
    try:
        supabase.table("quarantined_data").update({
            "status": "rejected",
            "reviewed_at": datetime.utcnow().isoformat()
        }).eq("id", item_id).execute()

        return {"status": "success", "message": "Item rejected"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/refinement/quarantined/{item_id}/ai-suggest")
async def generate_ai_suggestion(item_id: str):
    """AI를 사용해 수정 제안 생성"""
    try:
        import google.generativeai as genai

        # 1. Get quarantined item
        res = supabase.table("quarantined_data").select("*").eq("id", item_id).execute()
        if not res.data:
            raise HTTPException(status_code=404, detail="Item not found")

        item = res.data[0]
        original_data = item.get("original_data", {})
        error_type = item.get("error_type")

        # 2. Generate AI suggestion
        if not settings.GOOGLE_API_KEY:
            raise HTTPException(status_code=500, detail="API key not configured")

        genai.configure(api_key=settings.GOOGLE_API_KEY)
        model = genai.GenerativeModel('gemini-2.5-flash')

        prompt = f"""Fix the following data issue.

Error Type: {error_type}
Original Data: {json.dumps(original_data, ensure_ascii=False)}

Provide a corrected version of the data in JSON format.
Only output the corrected JSON, no explanation."""

        response = await asyncio.get_event_loop().run_in_executor(
            None, lambda: model.generate_content(prompt)
        )

        # Parse response
        try:
            suggested_fix = json.loads(response.text.strip())
        except:
            suggested_fix = {"ai_suggestion": response.text.strip()}

        # 3. Update quarantined item with suggestion
        supabase.table("quarantined_data").update({
            "suggested_fix": suggested_fix,
            "updated_at": datetime.utcnow().isoformat()
        }).eq("id", item_id).execute()

        return {"status": "success", "suggestion": suggested_fix}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/refinement/batch-fix")
async def run_batch_ai_fix(background_tasks: BackgroundTasks, type: str = "target_normalization"):
    """배치 AI 수정 작업 실행"""
    try:
        job_id = f"batch_{type}_{uuid4().hex[:8]}"

        # Create job record (if table exists)
        try:
            supabase.table("sync_jobs").insert({
                "id": job_id,
                "status": "queued",
                "source": f"batch_fix_{type}",
                "started_at": datetime.utcnow().isoformat()
            }).execute()
        except:
            pass

        # For now, return immediately (actual batch processing would be in background)
        return {
            "status": "started",
            "jobId": job_id,
            "type": type,
            "message": f"Batch {type} job started"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/refinement/health")
async def get_data_health():
    """데이터 건강 상태 조회"""
    try:
        # Use view if exists, otherwise calculate manually
        metrics = []

        tables = [
            ("antibody_library", "target_normalized", "embedding"),
            ("commercial_reagents", "target_normalized", "smiles"),
            ("golden_set", "target", "orr_pct")
        ]

        for table_name, norm_field, extra_field in tables:
            try:
                total_res = supabase.table(table_name).select("id", count="exact").execute()
                total = total_res.count or 0

                norm_res = supabase.table(table_name).select("id", count="exact").not_.is_(norm_field, "null").execute()
                normalized = norm_res.count or 0

                extra_res = supabase.table(table_name).select("id", count="exact").not_.is_(extra_field, "null").execute()
                extra_count = extra_res.count or 0

                metrics.append({
                    "sourceTable": table_name,
                    "totalCount": total,
                    "normalizedCount": normalized,
                    "embeddedCount": extra_count,
                    "normalizationRate": round((normalized / total * 100) if total > 0 else 0, 1),
                    "embeddingRate": round((extra_count / total * 100) if total > 0 else 0, 1)
                })
            except Exception as table_error:
                print(f"Error fetching health for {table_name}: {table_error}")
                continue

        return {"metrics": metrics}
    except Exception as e:
        print(f"Data health error: {e}")
        return {"metrics": []}


@router.get("/refinement/jobs")
async def get_batch_jobs(limit: int = 10):
    """배치 작업 목록 조회"""
    try:
        res = supabase.table("sync_jobs").select("*").order("started_at", desc=True).limit(limit).execute()

        jobs = []
        for item in (res.data or []):
            jobs.append({
                "id": item.get("id"),
                "type": item.get("source", "unknown").replace("batch_fix_", ""),
                "status": item.get("status"),
                "totalItems": item.get("total_items", 0),
                "processedItems": item.get("processed_items", 0),
                "successCount": item.get("success_count", 0),
                "errorCount": item.get("error_count", 0),
                "startedAt": item.get("started_at"),
                "completedAt": item.get("completed_at")
            })

        return {"jobs": jobs}
    except Exception as e:
        print(f"Jobs fetch error: {e}")
        return {"jobs": []}


# ============================================================
# v2.2 Management - Audit & Lineage
# ============================================================

@router.get("/audit/logs")
async def get_audit_logs(
    category: str = "all",
    action: str = "all",
    dateRange: str = "7d",
    limit: int = 50,
    offset: int = 0
):
    """관리자 감사 로그 조회"""
    try:
        query = supabase.table("admin_audit_logs").select("*", count="exact")

        if category != "all":
            query = query.eq("action_category", category)
        if action != "all":
            query = query.eq("action_type", action)

        # Date range filter
        from datetime import timedelta
        days = int(dateRange.replace("d", ""))
        cutoff = (datetime.utcnow() - timedelta(days=days)).isoformat()
        query = query.gte("created_at", cutoff)

        res = query.order("created_at", desc=True).range(offset, offset + limit - 1).execute()

        logs = []
        for item in (res.data or []):
            logs.append({
                "id": item.get("id"),
                "adminId": item.get("admin_id"),
                "adminEmail": item.get("admin_email", "admin@astraforge.io"),
                "actionType": item.get("action_type"),
                "actionCategory": item.get("action_category"),
                "targetType": item.get("target_type"),
                "targetId": item.get("target_id"),
                "beforeValue": item.get("before_value"),
                "afterValue": item.get("after_value"),
                "reason": item.get("reason"),
                "ipAddress": item.get("ip_address", "127.0.0.1"),
                "userAgent": item.get("user_agent", ""),
                "createdAt": item.get("created_at")
            })

        return {"logs": logs, "total": res.count or 0}
    except Exception as e:
        print(f"Audit logs error: {e}")
        return {"logs": [], "total": 0}


# ============================================================
# v2.2 Management - Analytics & Cost
# ============================================================

@router.get("/analytics")
async def get_analytics(range: str = "7d"):
    """분석 데이터 조회"""
    try:
        from datetime import timedelta

        days = int(range.replace("d", ""))
        cutoff = (datetime.utcnow() - timedelta(days=days)).isoformat()

        # Cost summary from budget_configs
        budget_res = supabase.table("budget_configs").select("*").execute()

        costs = []
        for item in (budget_res.data or []):
            costs.append({
                "provider": item.get("api_provider"),
                "dailyCost": float(item.get("current_usage_usd", 0)),
                "monthlyCost": float(item.get("current_monthly_usd", 0)),
                "dailyLimit": float(item.get("daily_limit_usd", 500)),
                "monthlyLimit": float(item.get("monthly_limit_usd", 10000)),
                "trend": 0  # Would need historical data to calculate
            })

        # Usage data from api_usage_logs (if exists)
        usage = []
        try:
            usage_res = supabase.table("api_usage_logs").select("*").gte("created_at", cutoff).order("created_at").execute()
            # Aggregate by day
            by_day = {}
            for log in (usage_res.data or []):
                day = log.get("created_at", "")[:10]
                if day not in by_day:
                    by_day[day] = {"totalCalls": 0, "totalTokens": 0, "totalCost": 0, "latencies": [], "errors": 0}
                by_day[day]["totalCalls"] += 1
                by_day[day]["totalTokens"] += log.get("tokens_input", 0) + log.get("tokens_output", 0)
                by_day[day]["totalCost"] += float(log.get("cost_usd", 0))
                by_day[day]["latencies"].append(log.get("latency_ms", 0))
                if not log.get("success"):
                    by_day[day]["errors"] += 1

            for day, data in by_day.items():
                avg_latency = sum(data["latencies"]) / len(data["latencies"]) if data["latencies"] else 0
                error_rate = (data["errors"] / data["totalCalls"] * 100) if data["totalCalls"] > 0 else 0
                usage.append({
                    "date": day,
                    "provider": "nvidia_nim",
                    "totalCalls": data["totalCalls"],
                    "totalTokens": data["totalTokens"],
                    "totalCost": data["totalCost"],
                    "avgLatency": int(avg_latency),
                    "errorRate": round(error_rate, 1)
                })
        except:
            pass

        # Session stats from design_runs (if exists)
        sessions = {
            "totalSessions": 0,
            "activeSessions": 0,
            "avgDuration": 0,
            "completionRate": 0,
            "byStage": {}
        }

        try:
            all_runs = supabase.table("design_runs").select("status, current_stage", count="exact").execute()
            sessions["totalSessions"] = all_runs.count or 0

            active_runs = supabase.table("design_runs").select("current_stage", count="exact").eq("status", "processing").execute()
            sessions["activeSessions"] = active_runs.count or 0

            # Group by stage
            for run in (active_runs.data or []):
                stage = run.get("current_stage", "unknown")
                sessions["byStage"][stage] = sessions["byStage"].get(stage, 0) + 1

            completed = supabase.table("design_runs").select("id", count="exact").eq("status", "completed").execute()
            if sessions["totalSessions"] > 0:
                sessions["completionRate"] = round((completed.count or 0) / sessions["totalSessions"] * 100, 1)
        except:
            pass

        return {
            "costs": costs,
            "usage": usage,
            "sessions": sessions
        }
    except Exception as e:
        print(f"Analytics error: {e}")
        return {"costs": [], "usage": [], "sessions": {}}



