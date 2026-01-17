"""
Clinical Agent - 임상 설계
임상 1상 프로토콜 초안 작성 및 권고사항 생성
"""
from typing import Dict, Any
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate

from app.agents.state import ADCState
from app.core.config import settings


CLINICAL_SYSTEM_PROMPT = """You are a Clinical Development Agent for ADC programs.

Your role:
1. Design Phase 1 dose-escalation protocol
2. Recommend starting dose based on preclinical data
3. Suggest biomarkers and endpoints
4. Identify key inclusion/exclusion criteria

Consider:
- Toxicity profile from previous analysis
- Target expression requirements
- Patient population optimization
- Companion diagnostic needs

Output format (JSON):
{
    "recommended_starting_dose": "1.2 mg/kg",
    "dose_escalation_scheme": "3+3 design",
    "primary_endpoints": ["MTD", "DLT", "RP2D"],
    "secondary_endpoints": ["ORR", "DCR", "DoR"],
    "biomarkers": ["Target expression (IHC)", "Circulating tumor DNA"],
    "patient_population": "Advanced solid tumors with target expression >10%",
    "key_exclusions": ["Prior ADC therapy", "Neuropathy Grade 2+"],
    "special_considerations": "<any notes>"
}
"""


async def run_clinical_agent(state: ADCState) -> Dict[str, Any]:
    """
    임상 기획 에이전트 실행
    
    Args:
        state: 현재 ADC 분석 상태
    
    Returns:
        Dict: 임상 프로토콜 권고사항
    """
    llm = ChatOpenAI(
        model="gpt-4o",
        temperature=0.2,
        api_key=settings.OPENAI_API_KEY
    )
    
    # 독성 정보 요약
    tox_summary = "Unknown risks"
    if state.toxicity_risks:
        tox_summary = ", ".join([f"{r.risk_type} ({r.severity})" for r in state.toxicity_risks[:3]])
    
    prompt = ChatPromptTemplate.from_messages([
        ("system", CLINICAL_SYSTEM_PROMPT),
        ("user", """Design Phase 1 protocol for the following ADC:

Target: {target_name}
Antibody: {antibody_type}
Payload: {payload_id}
DAR: {dar}

Known Toxicity Risks: {tox_summary}
Structure Stability: {stability}

Provide clinical development recommendations.""")
    ])
    
    chain = prompt | llm
    
    stability = "Unknown"
    if state.structure_analysis:
        stability = f"{state.structure_analysis.stability_score}/100, Aggregation: {state.structure_analysis.aggregation_risk}"
    
    response = await chain.ainvoke({
        "target_name": state.input.target_name,
        "antibody_type": state.input.antibody_type,
        "payload_id": state.input.payload_id,
        "dar": state.input.dar,
        "tox_summary": tox_summary,
        "stability": stability
    })
    
    # TODO: JSON 파싱
    # 현재는 Mock 데이터 반환
    protocol = {
        "recommended_starting_dose": "0.8 mg/kg" if state.input.payload_id == "mmae" else "1.2 mg/kg",
        "dose_escalation_scheme": "3+3 design",
        "primary_endpoints": ["MTD", "DLT", "RP2D"],
        "secondary_endpoints": ["ORR", "DCR", "DoR", "PFS"],
        "biomarkers": [
            f"{state.input.target_name} expression (IHC 2+/3+)",
            "Circulating tumor DNA",
            "Serum ADC concentration"
        ],
        "patient_population": f"Advanced solid tumors with {state.input.target_name} expression",
        "key_exclusions": [
            "Prior ADC therapy within 6 months",
            "Neuropathy Grade 2+" if state.input.payload_id == "mmae" else "ILD history",
            "LVEF < 50%"
        ],
        "estimated_patients": "30-50 patients",
        "estimated_duration": "12-18 months"
    }
    
    return protocol
