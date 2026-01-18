"""
Patent Agent - 특허 분석
Tavily/Google Patents 연동하여 IP 리스크 평가
"""
from typing import List
from langchain_core.prompts import ChatPromptTemplate

from app.agents.state import ADCState, PatentInfo
from app.core.gemini import get_gemini_model
from app.core.config import settings


PATENT_SYSTEM_PROMPT = """You are a Patent Analysis Agent for ADC development.

Your role:
1. Search for relevant patents related to the ADC configuration
2. Identify potential IP conflicts or FTO (Freedom-to-Operate) issues
3. Assess patent expiry dates and their implications
4. Recommend IP strategy

Key patent areas to consider:
- Antibody composition patents
- Linker chemistry patents
- Payload molecule patents
- Conjugation method patents
- Target-specific ADC patents

Output format (JSON array):
[
    {
        "patent_id": "US10123456",
        "title": "...",
        "status": "Active|Expired|Pending",
        "expiry_date": "2025-06-15",
        "risk_level": "Low|Medium|High"
    },
    ...
]
"""


async def run_patent_agent(state: ADCState) -> List[PatentInfo]:
    """
    특허 분석 에이전트 실행 (Gemini 2.0 Flash 사용)
    
    Args:
        state: 현재 ADC 분석 상태
    
    Returns:
        List[PatentInfo]: 관련 특허 목록
    """
    # Gemini 2.0 Flash 사용 (빠른 분석용)
    llm = get_gemini_model(temperature=0)
    
    prompt = ChatPromptTemplate.from_messages([
        ("system", PATENT_SYSTEM_PROMPT),
        ("user", """Analyze patent landscape for the following ADC:

Target: {target_name}
Antibody: {antibody_type}
Payload: {payload_id}
Linker: {linker_id}

Identify relevant patents and IP risks.""")
    ])
    
    chain = prompt | llm
    
    response = await chain.ainvoke({
        "target_name": state.input.target_name,
        "antibody_type": state.input.antibody_type,
        "payload_id": state.input.payload_id,
        "linker_id": state.input.linker_id
    })
    
    # TODO: Tavily API 연동으로 실제 특허 검색
    # 현재는 Mock 데이터 반환
    patents = [
        PatentInfo(
            patent_id="US9688761",
            title="Antibody-drug conjugates with improved stability",
            status="Expired",
            expiry_date="2024-03-15",
            risk_level="Low"
        ),
        PatentInfo(
            patent_id="US10570205",
            title=f"{state.input.payload_id.upper()} conjugation methods",
            status="Active",
            expiry_date="2032-08-22",
            risk_level="Medium" if state.input.payload_id == "mmae" else "Low"
        ),
    ]
    
    return patents
