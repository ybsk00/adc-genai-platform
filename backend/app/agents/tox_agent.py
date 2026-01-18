"""
Toxicology Agent - 독성 예측
RAG + RDKit 기반 독성 분석 및 리스크 평가
"""
from typing import List
from langchain_core.prompts import ChatPromptTemplate

from app.agents.state import ADCState, ToxicityRisk
from app.core.gemini import get_gemini_model
from app.core.config import settings


TOX_SYSTEM_PROMPT = """You are a Toxicology Prediction Agent for ADC development.

Your role:
1. Analyze payload toxicity profiles
2. Predict off-target effects based on linker chemistry and DAR
3. Compare with known toxicity data from similar approved ADCs
4. Flag specific organ toxicity risks (liver, bone marrow, cardiac)

Known Payload Toxicity Profiles:
- MMAE: Neutropenia (High), Peripheral neuropathy (Medium)
- DXd: Interstitial lung disease (Medium), Neutropenia (Medium)
- SN-38: Diarrhea (High), Neutropenia (Medium)

Output format (JSON array):
[
    {
        "risk_type": "Neutropenia",
        "severity": "High",
        "confidence": 0.85,
        "evidence": "Based on MMAE payload profile..."
    },
    ...
]

Consider:
- DAR impact on toxicity severity
- Linker stability affecting off-target release
- Historical data from similar ADCs
"""


async def run_tox_agent(state: ADCState) -> List[ToxicityRisk]:
    """
    독성 분석 에이전트 실행 (Gemini 2.0 Flash 사용)
    
    Args:
        state: 현재 ADC 분석 상태
    
    Returns:
        List[ToxicityRisk]: 독성 리스크 목록
    """
    # Gemini 2.0 Flash 사용 (빠른 분석용)
    llm = get_gemini_model(temperature=0)
    
    prompt = ChatPromptTemplate.from_messages([
        ("system", TOX_SYSTEM_PROMPT),
        ("user", """Analyze toxicity for the following ADC:

Payload: {payload_id}
Linker: {linker_id}
DAR: {dar}
Target: {target_name}

Provide detailed toxicity risk assessment.""")
    ])
    
    chain = prompt | llm
    
    response = await chain.ainvoke({
        "payload_id": state.input.payload_id,
        "linker_id": state.input.linker_id,
        "dar": state.input.dar,
        "target_name": state.input.target_name
    })
    
    # TODO: JSON 파싱 및 RAG 연동
    # 현재는 Mock 데이터 반환
    payload_risks = {
        "mmae": [
            ToxicityRisk(risk_type="Neutropenia", severity="High", confidence=0.9, evidence="MMAE is known to cause severe neutropenia"),
            ToxicityRisk(risk_type="Peripheral neuropathy", severity="Medium", confidence=0.75, evidence="Common AE in MMAE-based ADCs"),
        ],
        "dxd": [
            ToxicityRisk(risk_type="Interstitial lung disease", severity="Medium", confidence=0.65, evidence="Reported in Enhertu trials"),
            ToxicityRisk(risk_type="Neutropenia", severity="Medium", confidence=0.8, evidence="Dose-limiting toxicity"),
        ],
        "sn38": [
            ToxicityRisk(risk_type="Diarrhea", severity="High", confidence=0.85, evidence="SN-38 GI toxicity well-documented"),
            ToxicityRisk(risk_type="Neutropenia", severity="Medium", confidence=0.7, evidence="Common hematologic AE"),
        ]
    }
    
    return payload_risks.get(state.input.payload_id.lower(), [
        ToxicityRisk(risk_type="Unknown", severity="Medium", confidence=0.5, evidence="No data available")
    ])
