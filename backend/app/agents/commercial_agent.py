from typing import Dict, Any
from langchain_openai import ChatOpenAI
from app.core.config import settings
from app.agents.state import ADCState, CommercialInfo

COMMERCIAL_SYSTEM_PROMPT = """You are a Procurement Specialist for ADC research.
Your goal is to check the commercial availability of ADC components from the Ambeed database.

Role:
1. Identify if the specific Payload and Linker are available for purchase.
2. Prioritize "Linker-Payload Conjugates" (pre-assembled) as they save synthesis time.
3. Estimate the total cost for 10mg of material.
4. Assess feasibility:
   - High (Conjugate available in stock)
   - Medium (Individual components available)
   - Low (Custom synthesis required)

Output JSON format:
{
    "feasibility_score": 90,
    "total_estimated_cost": "$450",
    "availability_status": "Immediate",
    "recommendation": "Buy Conjugate Cat# A1234",
    "components": {
        "payload": {"name": "MMAE", "cat": "A100", "price": "$100/5mg"},
        "conjugate": {"name": "Mc-MMAE", "cat": "A200", "price": "$300/5mg"}
    }
}
"""

async def run_commercial_agent(state: ADCState) -> CommercialInfo:
    """
    상용화 분석 에이전트
    Ambeed DB를 검색하여 시약 구매 가능성 및 비용 견적을 산출합니다.
    비용 절감을 위해 FAST_LLM (gpt-4o-mini)을 사용합니다.
    """
    
    # [Cost Saving] 단순 검색이므로 gpt-4o-mini 사용
    llm = ChatOpenAI(
        model=settings.FAST_LLM,
        temperature=0,
        api_key=settings.OPENAI_API_KEY
    )
    
    # TODO: 실제 RAG 검색 로직 (Ambeed DB 조회) 구현 필요
    # 현재는 Mock Data로 동작
    
    payload_id = state.input.payload_id.lower()
    linker_id = state.input.linker_id.lower()
    
    # Mock Data: MMAE / DXd / SN-38
    if "mmae" in payload_id:
        return CommercialInfo(
            feasibility_score=95.0,
            total_estimated_cost="$231.00 (5mg Conjugate)",
            availability="In Stock",
            conjugate_info={
                "name": "Mc-VC-PAB-MMAE",
                "cat_no": "A503318",
                "price": "$231.00",
                "supplier": "Ambeed",
                "url": "https://www.ambeed.com/products/A503318.html"
            },
            payload_info={
                "name": "MMAE",
                "cat_no": "A100",
                "price": "$150.00",
                "supplier": "Ambeed"
            }
        )
    elif "dxd" in payload_id:
        return CommercialInfo(
            feasibility_score=85.0,
            total_estimated_cost="$450.00 (Custom Synthesis)",
            availability="Limited Stock",
            conjugate_info=None,
            payload_info={
                "name": "Deruxtecan (DXd)",
                "cat_no": "A200",
                "price": "$450.00",
                "supplier": "Ambeed"
            }
        )
    # [V2.0] Total Cost Estimation Logic (Commented out as requested)
    """
    payload_cost = float(payload_info.get("price", "0").replace("$", "").split("/")[0]) if payload_info else 0
    antibody_cost = 0
    target = state.input.target.lower()
    
    if "her2" in target:
        antibody_cost = 450  # Est. $450/mg (Trastuzumab Biosimilar)
    elif "trop" in target:
        antibody_cost = 500  # Est. $500/mg
    else:
        antibody_cost = 600  # Custom
        
    total_estimate = payload_cost + antibody_cost
    # total_estimated_cost = f"${total_estimate:.2/mg}"
    """

    return CommercialInfo(
        feasibility_score=40.0,
        total_estimated_cost="Unknown (Quote Required)",
        availability="Custom Synthesis Required",
        payload_info={"name": payload_id, "status": "Not Found"},
        linker_info={"name": linker_id, "status": "Not Found"}
    )
