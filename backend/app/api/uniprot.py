"""
UniProt API Integration
Target Antigen 검색 및 자동완성을 위한 UniProt API 연동
"""
from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel
from typing import List, Optional
import httpx
import logging
from functools import lru_cache
from datetime import datetime, timedelta

logger = logging.getLogger(__name__)

router = APIRouter()

# UniProt API Base URL
UNIPROT_API_BASE = "https://rest.uniprot.org/uniprotkb"

# 캐시 만료 시간 (1시간)
CACHE_TTL = timedelta(hours=1)
_cache: dict = {}
_cache_timestamps: dict = {}


class AntigenResult(BaseModel):
    """항원 검색 결과"""
    uniprot_id: str
    name: str
    gene_name: Optional[str] = None
    organism: str
    description: Optional[str] = None
    sequence_length: Optional[int] = None


class AntigenSearchResponse(BaseModel):
    """항원 검색 응답"""
    query: str
    results: List[AntigenResult]
    total_count: int


def _get_cached(key: str) -> Optional[dict]:
    """캐시에서 값 조회 (만료 확인)"""
    if key in _cache:
        if datetime.now() - _cache_timestamps[key] < CACHE_TTL:
            return _cache[key]
        else:
            del _cache[key]
            del _cache_timestamps[key]
    return None


def _set_cache(key: str, value: dict):
    """캐시에 값 저장"""
    _cache[key] = value
    _cache_timestamps[key] = datetime.now()


@router.get("/search", response_model=AntigenSearchResponse)
async def search_antigens(
    query: str = Query(..., min_length=2, description="검색어 (최소 2자)"),
    limit: int = Query(10, ge=1, le=50, description="결과 개수 제한"),
    organism: str = Query("human", description="생물종 필터 (human, mouse, all)")
):
    """
    UniProt에서 항원 검색

    Target Antigen 입력 시 자동완성을 위한 API

    - **query**: 검색어 (예: HER2, TROP2, EGFR)
    - **limit**: 최대 결과 개수 (기본 10개)
    - **organism**: 생물종 필터 (human, mouse, all)
    """
    # 캐시 키 생성
    cache_key = f"{query}:{limit}:{organism}"
    cached = _get_cached(cache_key)
    if cached:
        logger.info(f"[UniProt] Cache hit for query: {query}")
        return cached

    # 생물종 필터 쿼리 생성
    organism_filter = ""
    if organism == "human":
        organism_filter = " AND (organism_id:9606)"
    elif organism == "mouse":
        organism_filter = " AND (organism_id:10090)"

    # UniProt 검색 쿼리 구성
    # gene_name, protein_name에서 검색
    search_query = f"(gene:{query}* OR protein_name:{query}*){organism_filter}"

    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            response = await client.get(
                f"{UNIPROT_API_BASE}/search",
                params={
                    "query": search_query,
                    "format": "json",
                    "size": limit,
                    "fields": "accession,protein_name,gene_names,organism_name,length"
                }
            )

            if response.status_code != 200:
                logger.error(f"[UniProt] API error: {response.status_code}")
                raise HTTPException(502, "UniProt API error")

            data = response.json()

            results = []
            for entry in data.get("results", []):
                # 단백질 이름 추출
                protein_name = "Unknown"
                if "proteinDescription" in entry:
                    rec_name = entry["proteinDescription"].get("recommendedName")
                    if rec_name and "fullName" in rec_name:
                        protein_name = rec_name["fullName"].get("value", "Unknown")

                # 유전자 이름 추출
                gene_name = None
                if "genes" in entry and len(entry["genes"]) > 0:
                    gene_data = entry["genes"][0]
                    if "geneName" in gene_data:
                        gene_name = gene_data["geneName"].get("value")

                # 생물종 추출
                organism_name = entry.get("organism", {}).get("scientificName", "Unknown")

                results.append(AntigenResult(
                    uniprot_id=entry.get("primaryAccession", ""),
                    name=protein_name,
                    gene_name=gene_name,
                    organism=organism_name,
                    sequence_length=entry.get("sequence", {}).get("length")
                ))

            response_data = AntigenSearchResponse(
                query=query,
                results=results,
                total_count=len(results)
            )

            # 캐시에 저장
            _set_cache(cache_key, response_data.model_dump())

            logger.info(f"[UniProt] Found {len(results)} results for query: {query}")
            return response_data

    except httpx.TimeoutException:
        logger.error(f"[UniProt] Timeout for query: {query}")
        raise HTTPException(504, "UniProt API timeout")
    except Exception as e:
        logger.exception(f"[UniProt] Error: {e}")
        raise HTTPException(500, f"Search failed: {str(e)}")


@router.get("/entry/{uniprot_id}")
async def get_antigen_details(uniprot_id: str):
    """
    UniProt 항원 상세 정보 조회

    - **uniprot_id**: UniProt Accession ID (예: P04626)
    """
    cache_key = f"entry:{uniprot_id}"
    cached = _get_cached(cache_key)
    if cached:
        return cached

    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            response = await client.get(
                f"{UNIPROT_API_BASE}/{uniprot_id}",
                params={
                    "format": "json",
                    "fields": "accession,protein_name,gene_names,organism_name,length,sequence,cc_function,cc_subcellular_location"
                }
            )

            if response.status_code == 404:
                raise HTTPException(404, f"UniProt entry not found: {uniprot_id}")

            if response.status_code != 200:
                raise HTTPException(502, "UniProt API error")

            entry = response.json()

            # 상세 정보 파싱
            protein_name = "Unknown"
            if "proteinDescription" in entry:
                rec_name = entry["proteinDescription"].get("recommendedName")
                if rec_name and "fullName" in rec_name:
                    protein_name = rec_name["fullName"].get("value", "Unknown")

            gene_name = None
            if "genes" in entry and len(entry["genes"]) > 0:
                gene_data = entry["genes"][0]
                if "geneName" in gene_data:
                    gene_name = gene_data["geneName"].get("value")

            # 기능 설명 추출
            function_text = None
            if "comments" in entry:
                for comment in entry["comments"]:
                    if comment.get("commentType") == "FUNCTION":
                        texts = comment.get("texts", [])
                        if texts:
                            function_text = texts[0].get("value")
                        break

            result = {
                "uniprot_id": entry.get("primaryAccession"),
                "name": protein_name,
                "gene_name": gene_name,
                "organism": entry.get("organism", {}).get("scientificName"),
                "sequence_length": entry.get("sequence", {}).get("length"),
                "function": function_text,
                "sequence": entry.get("sequence", {}).get("value")
            }

            _set_cache(cache_key, result)
            return result

    except HTTPException:
        raise
    except Exception as e:
        logger.exception(f"[UniProt] Error fetching entry: {e}")
        raise HTTPException(500, f"Failed to fetch entry: {str(e)}")


# 자주 사용되는 ADC 타겟 항원 목록 (오프라인 fallback)
COMMON_ADC_TARGETS = [
    {"gene_name": "HER2", "name": "Receptor tyrosine-protein kinase erbB-2", "uniprot_id": "P04626", "organism": "Homo sapiens"},
    {"gene_name": "TROP2", "name": "Tumor-associated calcium signal transducer 2", "uniprot_id": "P09758", "organism": "Homo sapiens"},
    {"gene_name": "EGFR", "name": "Epidermal growth factor receptor", "uniprot_id": "P00533", "organism": "Homo sapiens"},
    {"gene_name": "CD33", "name": "Myeloid cell surface antigen CD33", "uniprot_id": "P20138", "organism": "Homo sapiens"},
    {"gene_name": "CD30", "name": "Tumor necrosis factor receptor superfamily member 8", "uniprot_id": "P28908", "organism": "Homo sapiens"},
    {"gene_name": "CD22", "name": "B-cell receptor CD22", "uniprot_id": "P20273", "organism": "Homo sapiens"},
    {"gene_name": "CD19", "name": "B-lymphocyte antigen CD19", "uniprot_id": "P15391", "organism": "Homo sapiens"},
    {"gene_name": "BCMA", "name": "Tumor necrosis factor receptor superfamily member 17", "uniprot_id": "Q02223", "organism": "Homo sapiens"},
    {"gene_name": "NECTIN4", "name": "Nectin-4", "uniprot_id": "Q96NY8", "organism": "Homo sapiens"},
    {"gene_name": "FRα", "name": "Folate receptor alpha", "uniprot_id": "P15328", "organism": "Homo sapiens"},
    {"gene_name": "PSMA", "name": "Glutamate carboxypeptidase 2", "uniprot_id": "Q04609", "organism": "Homo sapiens"},
    {"gene_name": "DLL3", "name": "Delta-like protein 3", "uniprot_id": "Q9NYJ7", "organism": "Homo sapiens"},
    {"gene_name": "CLDN18.2", "name": "Claudin-18", "uniprot_id": "P56856", "organism": "Homo sapiens"},
    {"gene_name": "GPNMB", "name": "Transmembrane glycoprotein NMB", "uniprot_id": "Q14956", "organism": "Homo sapiens"},
    {"gene_name": "MSLN", "name": "Mesothelin", "uniprot_id": "Q13421", "organism": "Homo sapiens"},
]


@router.get("/common-targets")
async def get_common_targets():
    """
    자주 사용되는 ADC 타겟 항원 목록

    오프라인 fallback 및 빠른 선택을 위한 사전 정의 목록
    """
    return {
        "targets": COMMON_ADC_TARGETS,
        "count": len(COMMON_ADC_TARGETS)
    }
