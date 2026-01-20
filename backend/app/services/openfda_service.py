"""
OpenFDA Service - 승인된 ADC 약물 라벨 데이터 수집
Payload/Linker 키워드로 역추적하여 골든셋의 'Approved' 데이터를 완성함
"""
import httpx
from typing import List, Dict, Any
import asyncio
from tenacity import retry, stop_after_attempt, wait_exponential

class OpenFDAService:
    BASE_URL = "https://api.fda.gov/drug/label.json"
    
    # [핵심] 3대 ADC 및 유명 독립군을 잡기 위한 '그물망' 키워드
    SEARCH_QUERIES = [
        # 1. Payload Suffix로 검색 (가장 강력함)
        'openfda.generic_name:"*vedotin"',      # MMAE/MMAF 계열 (Adcetris, Padcev, Polivy...)
        'openfda.generic_name:"*deruxtecan"',   # DXd 계열 (Enhertu)
        'openfda.generic_name:"*govitecan"',    # SN-38 계열 (Trodelvy)
        'openfda.generic_name:"*tansine"',      # DM1 계열 (Kadcyla)
        'openfda.generic_name:"*ozogamicin"',   # Calicheamicin 계열 (Mylotarg)
        'openfda.generic_name:"*soravtansine"', # DM4 계열 (Elahere)
        
        # 2. 약효 분류로 검색 (안전망)
        'openfda.pharm_class_epc:"Antibody-Drug Conjugate"'
    ]

    @retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
    async def fetch_labels(self, query: str, limit: int = 10) -> List[Dict[str, Any]]:
        """OpenFDA 라벨 데이터 검색"""
        params = {
            "search": query,
            "limit": limit
        }
        async with httpx.AsyncClient() as client:
            response = await client.get(self.BASE_URL, params=params, timeout=30.0)
            if response.status_code == 404:
                return [] # 검색 결과 없음
            response.raise_for_status()
            return response.json().get("results", [])

    def extract_golden_info(self, label: Dict[str, Any]) -> Dict[str, Any]:
        """
        라벨 데이터에서 '골든셋'에 필요한 핵심 정보만 정제
        """
        openfda = label.get("openfda", {})
        
        # 1. 기본 정보
        brand_name = openfda.get("brand_name", ["Unknown"])[0]
        generic_name = openfda.get("generic_name", ["Unknown"])[0]
        manufacturer = openfda.get("manufacturer_name", ["Unknown"])[0]
        
        # 2. 독성 정답지 (Boxed Warning)
        boxed_warning = label.get("boxed_warning", [])
        toxicity_summary = " ".join(boxed_warning) if boxed_warning else "No Boxed Warning"
        
        # 3. 구조 정보 (Description 섹션에서 링커/페이로드 텍스트 추출 시도)
        description = label.get("description", [])
        desc_text = " ".join(description) if description else ""
        
        return {
            "name": brand_name,
            "generic_name": generic_name,
            "category": "approved_drug",
            "manufacturer": manufacturer,
            "description": f"[{manufacturer}] {generic_name}. {toxicity_summary[:200]}...",
            "properties": {
                "generic_name": generic_name,
                "target": "Unknown (Check Description)", # NLP나 LLM으로 추가 추출 필요
                "payload_type": "Unknown (Check Generic Name)",
                "toxicity_profile": toxicity_summary, # 이게 진짜 독성 정답지
                "approval_status": "Approved",
                "full_description": desc_text
            },
            "status": "approved", # FDA 데이터는 무조건 승인된 것
            "outcome_type": "Success",
            "enrichment_source": "openfda"
        }

openfda_service = OpenFDAService()
