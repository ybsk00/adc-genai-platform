from fastapi import APIRouter
from typing import Optional
from app.core.supabase import supabase

router = APIRouter()

@router.get("/stats")
async def get_library_stats():
    """데이터 소스별 통계 조회"""
    stats = {}
    
    try:
        # 1. Clinical Trials (from golden_set_library where source='clinical_trials')
        clinical = supabase.table("golden_set_library").select("count", count="exact").eq("enrichment_source", "clinical_trials").execute()
        stats["clinical"] = clinical.count
        
        # 2. PubMed (from knowledge_base where source_type='PubMed')
        pubmed = supabase.table("knowledge_base").select("count", count="exact").eq("source_type", "PubMed").execute()
        stats["pubmed"] = pubmed.count
        
        # 3. Golden Set (Total in golden_set_library)
        goldenset = supabase.table("golden_set_library").select("count", count="exact").execute()
        stats["goldenset"] = goldenset.count
        
        # 4. News (Placeholder)
        stats["news"] = 0
        
    except Exception as e:
        print(f"Stats Error: {e}")
        
    return stats

@router.get("/goldenset")
async def search_golden_set(
    target: Optional[str] = None,
    category: Optional[str] = None,
    review_required: Optional[bool] = None,
    page: int = 1,
    limit: int = 20
):
    """골든셋 목록 검색 (Real DB)"""
    try:
        query = supabase.table("golden_set_library").select("*", count="exact")
        
        if target:
            query = query.ilike("properties->>target", f"%{target}%")
            
        if review_required is not None:
            query = query.eq("review_required", review_required)
            
        # Pagination
            
        # Pagination
        start = (page - 1) * limit
        end = start + limit - 1
        query = query.range(start, end)
        
        result = query.execute()
        
        return {
            "items": result.data,
            "total": result.count,
            "page": page,
            "limit": limit
        }
    except Exception as e:
        print(f"Search Error: {e}")
        return {"items": [], "total": 0, "page": page, "limit": limit}
