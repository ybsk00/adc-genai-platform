"""
Structure Agent - 3D 구조 분석
AlphaFold/BioNeMo 연동하여 항체 구조 분석 및 접합 부위 식별
"""
from typing import Optional
from langchain_core.prompts import ChatPromptTemplate

from app.agents.state import ADCState, StructureAnalysis
from app.core.gemini import get_gemini_model
from app.core.config import settings


STRUCTURE_SYSTEM_PROMPT = """You are a Structure Analysis Agent specializing in ADC (Antibody-Drug Conjugate) molecular structures.

Your role:
1. Analyze 3D protein structures from the given antibody information
2. Identify optimal conjugation sites on the antibody
3. Predict structural stability after payload attachment
4. Assess aggregation risk based on DAR and linker chemistry

Output format (JSON):
{
    "stability_score": <0-100>,
    "conjugation_sites": ["Cys-239", "Cys-242", ...],
    "aggregation_risk": "Low|Medium|High",
    "analysis_notes": "<detailed analysis>"
}

Consider:
- DAR (Drug-to-Antibody Ratio) impact on aggregation
- Linker type effect on conjugation efficiency
- Known antibody variant stability issues
"""


async def run_structure_agent(state: ADCState) -> StructureAnalysis:
    """
    구조 분석 에이전트 실행 (Gemini 2.0 Flash 사용)
    
    Args:
        state: 현재 ADC 분석 상태
    
    Returns:
        StructureAnalysis: 구조 분석 결과
    """
    # Gemini 2.0 Flash 사용 (빠른 분석용)
    llm = get_gemini_model(temperature=0)
    
    prompt = ChatPromptTemplate.from_messages([
        ("system", STRUCTURE_SYSTEM_PROMPT),
        ("user", """Analyze the following ADC configuration:

Target: {target_name}
Antibody: {antibody_type}
Custom Sequence: {custom_sequence}
Payload: {payload_id}
Linker: {linker_id}
DAR: {dar}

Provide a comprehensive structural analysis.""")
    ])
    
    chain = prompt | llm
    
    response = await chain.ainvoke({
        "target_name": state.input.target_name,
        "antibody_type": state.input.antibody_type,
        "custom_sequence": state.input.custom_sequence or "N/A",
        "payload_id": state.input.payload_id,
        "linker_id": state.input.linker_id,
        "dar": state.input.dar
    })
    
    # [Mock Logic] 유사도 계산 시뮬레이션
    is_similar_to_enhertu = (
        "her2" in state.input.target_name.lower() and 
        "dxd" in state.input.payload_id.lower()
    )
    
    comparison_note = ""
    if is_similar_to_enhertu:
        comparison_note = "Structure aligns with Trastuzumab deruxtecan (RMSD < 0.5Å). Conjugation sites match Enhertu profile."
    else:
        comparison_note = "Standard IgG1 scaffold structure. No significant deviation from baseline."

    return StructureAnalysis(
        stability_score=92.5 if is_similar_to_enhertu else 75.0,
        conjugation_sites=["Cys-239", "Cys-242", "Cys-265"],
        aggregation_risk="Low" if state.input.dar <= 4 else "Medium",
        pdb_data=None,  # AlphaFold 연동 시 채움
        analysis_notes=comparison_note
    )

