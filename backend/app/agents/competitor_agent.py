"""
Competitor Agent - 경쟁사 분석
Perplexity API 연동하여 시장 동향 및 경쟁사 현황 분석
"""
from typing import List
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate

from app.agents.state import ADCState, CompetitorInfo
from app.core.config import settings


COMPETITOR_SYSTEM_PROMPT = """You are a Competitive Intelligence Agent for ADC development.

Your role:
1. Identify competitors developing similar ADCs
2. Track development stages (Phase 1, 2, 3, Approved)
3. Analyze market positioning and differentiation opportunities
4. Provide strategic insights

Key information to gather:
- Active clinical trials with similar targets
- Recent regulatory approvals
- Pipeline announcements
- M&A activities in the ADC space

Output format (JSON array):
[
    {
        "company": "Daiichi Sankyo",
        "drug_name": "Enhertu",
        "phase": "Approved",
        "target": "HER2",
        "notes": "First-in-class for HER2-low"
    },
    ...
]
"""


async def run_competitor_agent(state: ADCState) -> List[CompetitorInfo]:
    """
    경쟁사 분석 에이전트 실행
    
    Args:
        state: 현재 ADC 분석 상태
    
    Returns:
        List[CompetitorInfo]: 경쟁사 정보 목록
    """
    llm = ChatOpenAI(
        model="gpt-4o",
        temperature=0,
        api_key=settings.OPENAI_API_KEY
    )
    
    prompt = ChatPromptTemplate.from_messages([
        ("system", COMPETITOR_SYSTEM_PROMPT),
        ("user", """Analyze competitive landscape for the following ADC:

Target: {target_name}
Payload type: {payload_id}

Identify active competitors and their development status.""")
    ])
    
    chain = prompt | llm
    
    response = await chain.ainvoke({
        "target_name": state.input.target_name,
        "payload_id": state.input.payload_id
    })
    
    # TODO: Perplexity API 연동
    # 현재는 Mock 데이터 반환
    target_competitors = {
        "her2": [
            CompetitorInfo(company="Daiichi Sankyo/AstraZeneca", drug_name="Enhertu", phase="Approved", target="HER2", notes="First-in-class for HER2-low"),
            CompetitorInfo(company="Roche", drug_name="Kadcyla", phase="Approved", target="HER2", notes="mAb + DM1"),
            CompetitorInfo(company="Seagen/Merck", drug_name="Tukysa+", phase="Phase 3", target="HER2", notes="Combination therapy"),
        ],
        "trop-2": [
            CompetitorInfo(company="Gilead", drug_name="Trodelvy", phase="Approved", target="TROP-2", notes="First TROP-2 ADC"),
            CompetitorInfo(company="Daiichi Sankyo", drug_name="DS-1062", phase="Phase 3", target="TROP-2", notes="DXd payload"),
        ],
        "liv-1": [
            CompetitorInfo(company="Seagen", drug_name="Ladiratuzumab vedotin", phase="Phase 2", target="LIV-1", notes="MMAE-based"),
        ]
    }
    
    target_lower = state.input.target_name.lower().replace(" ", "").replace("-", "")
    
    for key, competitors in target_competitors.items():
        if key in target_lower:
            return competitors
    
    # 기본 경쟁사
    return [
        CompetitorInfo(company="Multiple", drug_name="Various", phase="Phase 1-3", target=state.input.target_name, notes="Emerging target with limited competition")
    ]
