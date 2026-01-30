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
@router.get("/golden-set")
async def search_golden_set(
    target: Optional[str] = None,
    category: Optional[str] = None,
    review_required: Optional[bool] = None,
    status: Optional[str] = None,
    search: Optional[str] = None,
    page: int = 1,
    limit: int = 20
):
    """
    골든셋 목록 검색 (Real DB)

    - target: 타겟 항원 필터 (target_1, target_2)
    - category: 카테고리 필터
    - status: 상태 필터 (draft, approved, rejected)
    - search: 이름 검색
    """
    try:
        query = supabase.table("golden_set_library").select(
            "id, name, category, status, target_1, target_2, outcome_type, "
            "dar, molecular_weight, canonical_smiles, payload_smiles, linker_smiles, "
            "linker_type, antibody_format, orr_pct, os_months, pfs_months, "
            "confidence_score, ai_refined, created_at, updated_at",
            count="exact"
        )

        # 검색 필터
        if search:
            query = query.or_(f"name.ilike.%{search}%,category.ilike.%{search}%")

        if target:
            # target_1 또는 target_2에서 검색
            query = query.or_(f"target_1.ilike.%{target}%,target_2.ilike.%{target}%")

        if category:
            query = query.eq("category", category)

        if status:
            query = query.eq("status", status)

        if review_required is not None:
            query = query.eq("review_required", review_required)

        # 정렬 및 페이지네이션
        start = (page - 1) * limit
        end = start + limit - 1
        query = query.order("updated_at", desc=True).range(start, end)

        result = query.execute()

        # 응답 형식 통일 (items 키 사용)
        return {
            "items": result.data,
            "total": result.count or 0,
            "page": page,
            "limit": limit
        }
    except Exception as e:
        print(f"Golden Set Search Error: {e}")
        return {"items": [], "total": 0, "page": page, "limit": limit}


@router.get("/antibodies")
async def get_antibodies(
    page: int = 1,
    limit: int = 20,
    search: Optional[str] = None,
    source: Optional[str] = None
):
    """항체 라이브러리 조회 (cat_no 검색 포함)"""
    try:
        query = supabase.table("antibody_library").select("*", count="exact")

        if search:
            # OR search on product_name, cat_no, related_disease
            query = query.or_(f"product_name.ilike.%{search}%,cat_no.ilike.%{search}%,related_disease.ilike.%{search}%")

        if source:
            query = query.eq("source_name", source)

        start = (page - 1) * limit
        end = start + limit - 1
        query = query.order("crawled_at", desc=True).range(start, end)

        result = query.execute()

        return {
            "data": result.data,
            "total": result.count,
            "page": page,
            "limit": limit
        }
    except Exception as e:
        print(f"Antibody Search Error: {e}")
        return {"data": [], "total": 0, "page": page, "limit": limit}


@router.get("/reagents")
async def get_reagents(
    page: int = 1,
    limit: int = 20,
    search: Optional[str] = None,
    category: Optional[str] = None,
    missing_smiles: Optional[bool] = None,
    ai_refined: Optional[bool] = None,
    manual_override: Optional[bool] = None,
    source: Optional[str] = None
):
    """
    시약 라이브러리 조회 (Status 필터링 포함)

    - category: 카테고리 필터 (payload, linker, antibody, target 등)
    - missing_smiles: SMILES가 없는 항목만
    - ai_refined: AI 정제된 항목만
    - manual_override: 수동 수정된 항목만
    - source: 소스 필터 (Ambeed, MedChem Express 등)
    """
    try:
        query = supabase.table("commercial_reagents").select(
            "id, ambeed_cat_no, product_name, category, cas_number, "
            "smiles_code, molecular_weight, formula, stock_status, "
            "target, source_name, ai_refined, is_manual_override, "
            "payload_smiles, linker_smiles, full_smiles, crawled_at",
            count="exact"
        )

        if search:
            # OR search on product_name, cas_number, ambeed_cat_no
            query = query.or_(f"product_name.ilike.%{search}%,cas_number.ilike.%{search}%,ambeed_cat_no.ilike.%{search}%")

        # Category Filter (중요: payload, linker 구분)
        if category:
            query = query.ilike("category", f"%{category}%")

        # Status Filters
        if missing_smiles:
            query = query.is_("smiles_code", "null")

        if ai_refined is not None:
            query = query.eq("ai_refined", ai_refined)

        if manual_override is not None:
            query = query.eq("is_manual_override", manual_override)

        if source:
            query = query.eq("source_name", source)

        start = (page - 1) * limit
        end = start + limit - 1
        query = query.order("crawled_at", desc=True).range(start, end)

        result = query.execute()

        # 응답 형식 통일: items 키 사용 (프론트엔드 호환)
        return {
            "items": result.data,
            "total": result.count or 0,
            "page": page,
            "limit": limit
        }
    except Exception as e:
        print(f"Reagent Search Error: {e}")
        return {"items": [], "total": 0, "page": page, "limit": limit}


@router.post("/reagents/{id}/autofill-smiles")
async def autofill_reagent_smiles(id: str):
    """
    PubChem에서 SMILES 자동 채우기 (단일 레코드)
    CAS 또는 Name으로 조회 후 Desalt + MW 검증
    """
    try:
        from app.services.pubchem_service import pubchem_service

        # 1. 레코드 조회
        res = supabase.table("commercial_reagents").select("*").eq("id", id).execute()
        if not res.data:
            return {"error": "Record not found"}

        record = res.data[0]

        # 2. PubChem 조회
        result = await pubchem_service.autofill_smiles(record)

        if "error" in result:
            return result

        # 3. DB 업데이트 (선택적 - 자동 저장 시)
        # 여기서는 결과만 반환하고, 프론트엔드에서 확인 후 저장
        return {
            "status": "success",
            "record_id": id,
            "current_smiles": record.get("smiles_code"),
            "suggested_smiles": result.get("smiles_code"),
            "desalted_smiles": result.get("desalted_smiles"),
            "molecular_weight": result.get("molecular_weight"),
            "calculated_mw": result.get("calculated_mw"),
            "mw_difference": result.get("mw_difference"),
            "validation_status": result.get("validation_status"),
            "pubchem_cid": result.get("pubchem_cid"),
            "formula": result.get("formula")
        }
    except Exception as e:
        print(f"Autofill SMILES Error: {e}")
        return {"error": str(e)}
