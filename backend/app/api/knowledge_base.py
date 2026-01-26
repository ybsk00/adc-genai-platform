"""
Knowledge Base API - 비정형 데이터 관리
뉴스, 논문 등 텍스트 데이터를 조회하고 RAG 인덱싱을 수행
"""
from fastapi import APIRouter, HTTPException, Depends, BackgroundTasks
from pydantic import BaseModel, Field
from typing import Optional, List
from datetime import datetime
from uuid import uuid4

from app.core.supabase import supabase
from app.core.config import settings
# from app.services.rag_service import index_document # TODO: RAG 서비스 구현 필요

router = APIRouter()

class KnowledgeBaseItem(BaseModel):
    id: int
    source_type: str
    title: str
    summary: Optional[str]
    relevance_score: float
    source_tier: int
    rag_status: str
    created_at: str

class IndexRequest(BaseModel):
    """지식 베이스 인덱싱 요청"""
    source_type: str = Field(..., description="News, PubMed, Report")
    title: str
    content: str
    summary: Optional[str] = None
    url: Optional[str] = None

@router.get("/", response_model=List[KnowledgeBaseItem])
async def get_knowledge_base_items(
    status: Optional[str] = None,
    limit: int = 50
):
    """
    지식 베이스 목록 조회
    status: pending, indexed, excluded
    """
    try:
        query = supabase.table("knowledge_base").select("*").order("created_at", desc=True).limit(limit)
        
        if status:
            query = query.eq("rag_status", status)
            
        response = query.execute()
        return response.data
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/index")
async def add_to_knowledge_base(req: IndexRequest, background_tasks: BackgroundTasks):
    """
    새로운 지식 추가 및 인덱싱 요청
    """
    try:
        # 1. DB 저장
        new_item = {
            "source_type": req.source_type,
            "title": req.title,
            "content": req.content,
            "summary": req.summary or req.content[:200],
            "relevance_score": 0.0, # 초기값
            "source_tier": 2, # 기본값
            "rag_status": "pending"
        }
        
        res = supabase.table("knowledge_base").insert(new_item).execute()
        
        if not res.data:
            raise HTTPException(status_code=500, detail="Failed to insert data")
            
        item_id = res.data[0]['id']
        
        # 2. 백그라운드 인덱싱 (TODO: 실제 구현)
        # background_tasks.add_task(index_document, item_id, req.content)
        
        return {
            "status": "success",
            "id": item_id,
            "message": "Data added to Knowledge Base. Indexing pending."
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

class KnowledgeItemUpdate(BaseModel):
    title: Optional[str] = None
    content: Optional[str] = None
    summary: Optional[str] = None
    ai_reasoning: Optional[str] = None
    relevance_score: Optional[float] = None
    properties: Optional[dict] = None
    rag_status: Optional[str] = None

@router.patch("/{item_id}")
async def patch_knowledge_item(item_id: int, updates: KnowledgeItemUpdate):
    """
    지식 베이스 항목 수정
    """
    try:
        # Filter out None values
        update_data = {k: v for k, v in updates.dict().items() if v is not None}
        
        if not update_data:
            return {"message": "No updates provided"}

        res = supabase.table("knowledge_base").update(update_data).eq("id", item_id).execute()
        
        if not res.data:
            raise HTTPException(status_code=404, detail="Item not found")
            
        return res.data[0]
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))
