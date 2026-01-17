from fastapi import APIRouter
from typing import Optional

router = APIRouter()


@router.get("/goldenset")
async def search_golden_set(
    target: Optional[str] = None,
    category: Optional[str] = None,
    page: int = 1,
    limit: int = 20
):
    """골든셋 목록 검색"""
    # TODO: DB 검색 및 벡터 유사도 검색
    return {
        "items": [
            {
                "id": "gs_1",
                "name": "MMAE",
                "category": "payload",
                "properties": {"logP": 3.2, "mw": 718}
            },
            {
                "id": "gs_2", 
                "name": "Trastuzumab",
                "category": "antibody",
                "properties": {"target": "HER2"}
            }
        ],
        "total": 100,
        "page": page,
        "limit": limit
    }
