"""
Perplexity Service - 데이터 보강 (Enrichment)
외부 데이터 수집 시 비어있는 필드를 Perplexity API로 채움
"""
import httpx
from typing import Optional, Dict, Any
import json

from app.core.config import settings


PERPLEXITY_API_URL = "https://api.perplexity.ai/chat/completions"


async def enrich_drug_data(raw_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Perplexity API를 활용하여 ADC 약물 데이터 보강
    
    비용 절감을 위해 필수 필드(linker, payload)가 없을 때만 호출
    
    Args:
        raw_data: 원본 데이터 (drug_name, target 등)
    
    Returns:
        보강된 데이터 + enrichment_source
    """
    drug_name = raw_data.get("drug_name", "")
    target = raw_data.get("target", "")
    
    # 필수 필드가 있으면 보강 불필요
    if raw_data.get("linker") and raw_data.get("payload"):
        return {
            **raw_data,
            "enrichment_source": "original",
            "enriched": False
        }
    
    # Perplexity에게 부족한 정보 요청
    prompt = f"""Find detailed information about the ADC (Antibody-Drug Conjugate) drug '{drug_name}' targeting '{target}'.

Return ONLY a valid JSON object with these fields:
{{
    "linker_type": "linker name and type (cleavable/non-cleavable)",
    "payload_type": "payload name (e.g., MMAE, DXd, SN-38)",
    "dar": "Drug-to-Antibody Ratio (number)",
    "company": "developing company name",
    "mechanism": "brief mechanism of action",
    "clinical_phase": "current development phase"
}}

If information is not available, use null for that field."""

    try:
        async with httpx.AsyncClient() as client:
            response = await client.post(
                PERPLEXITY_API_URL,
                headers={
                    "Authorization": f"Bearer {settings.PERPLEXITY_API_KEY}",
                    "Content-Type": "application/json"
                },
                json={
                    "model": "sonar-medium-chat",
                    "messages": [
                        {"role": "user", "content": prompt}
                    ],
                    "temperature": 0.1,
                    "max_tokens": 500
                },
                timeout=30.0
            )
            response.raise_for_status()
            
            result = response.json()
            content = result.get("choices", [{}])[0].get("message", {}).get("content", "")
            
            # JSON 파싱 시도
            try:
                # 마크다운 코드 블록 제거
                if "```json" in content:
                    content = content.split("```json")[1].split("```")[0]
                elif "```" in content:
                    content = content.split("```")[1].split("```")[0]
                
                enriched_info = json.loads(content.strip())
                
                # 원본 데이터와 병합
                enriched_data = {
                    **raw_data,
                    "linker": enriched_info.get("linker_type") or raw_data.get("linker"),
                    "payload": enriched_info.get("payload_type") or raw_data.get("payload"),
                    "dar": enriched_info.get("dar") or raw_data.get("dar"),
                    "company": enriched_info.get("company") or raw_data.get("company"),
                    "enrichment_source": "perplexity_sonar_medium",
                    "enriched": True,
                    "enrichment_details": enriched_info
                }
                
                return enriched_data
                
            except json.JSONDecodeError:
                # JSON 파싱 실패 시 원본 반환
                return {
                    **raw_data,
                    "enrichment_source": "perplexity_failed",
                    "enriched": False,
                    "enrichment_error": "Failed to parse Perplexity response"
                }
                
    except httpx.HTTPError as e:
        return {
            **raw_data,
            "enrichment_source": "perplexity_error",
            "enriched": False,
            "enrichment_error": str(e)
        }


async def search_competitor_info(target: str, payload_type: str) -> Optional[Dict[str, Any]]:
    """
    특정 타겟/페이로드 조합의 경쟁사 정보 검색
    """
    prompt = f"""List all ADC drugs targeting '{target}' or using '{payload_type}' payload.

For each drug, provide:
- Drug name
- Company
- Current development phase
- Approval status (if approved)

Return as a JSON array."""

    try:
        async with httpx.AsyncClient() as client:
            response = await client.post(
                PERPLEXITY_API_URL,
                headers={
                    "Authorization": f"Bearer {settings.PERPLEXITY_API_KEY}",
                    "Content-Type": "application/json"
                },
                json={
                    "model": "sonar-medium-chat",
                    "messages": [{"role": "user", "content": prompt}],
                    "temperature": 0.1
                },
                timeout=30.0
            )
            response.raise_for_status()
            
            result = response.json()
            content = result.get("choices", [{}])[0].get("message", {}).get("content", "")
            
            # JSON 파싱
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0]
            elif "```" in content:
                content = content.split("```")[1].split("```")[0]
            
            return json.loads(content.strip())
            
    except Exception:
        return None
