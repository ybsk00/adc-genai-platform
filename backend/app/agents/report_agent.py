"""
Report Agent - 최종 리포트 생성
모든 에이전트 결과를 종합하여 Executive Summary 및 PDF 생성
"""
from typing import Dict, Any
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate

from app.agents.state import ADCState
from app.core.config import settings


REPORT_SYSTEM_PROMPT = """You are a Report Writer for ADC analysis results.

Your role:
1. Synthesize findings from all analysis domains
2. Generate an executive summary in Korean
3. Assign an overall grade (A, B+, B, C, D)
4. Provide a clear recommendation (Go, Conditional Go, No Go)

Grading criteria:
- A: Excellent profile, minimal risks, strong commercial potential
- B+: Good profile with manageable risks
- B: Acceptable profile with notable risks requiring attention
- C: Significant concerns, proceed with caution
- D: Major issues, not recommended to proceed

Output format (JSON):
{
    "grade": "B+",
    "recommendation": "Conditional Go",
    "summary": "<Korean executive summary, 3-5 sentences>",
    "key_strengths": ["...", "..."],
    "key_risks": ["...", "..."],
    "next_steps": ["...", "..."],
    "scores": {
        "efficacy": <0-100>,
        "toxicity": <0-100, higher is safer>,
        "properties": <0-100>,
        "patent": <0-100>,
        "market": <0-100>
    }
}
"""


async def run_report_agent(state: ADCState) -> Dict[str, Any]:
    """
    리포트 생성 에이전트 실행
    
    Args:
        state: 현재 ADC 분석 상태 (모든 에이전트 결과 포함)
    
    Returns:
        Dict: 최종 리포트 데이터
    """
    llm = ChatOpenAI(
        model="gpt-4o",
        temperature=0.3,
        api_key=settings.OPENAI_API_KEY
    )
    
    # 각 분석 결과 요약
    structure_summary = "분석 실패"
    if state.structure_analysis:
        structure_summary = f"안정성 {state.structure_analysis.stability_score}/100, 응집 리스크: {state.structure_analysis.aggregation_risk}"
    
    tox_summary = "분석 실패"
    if state.toxicity_risks:
        high_risks = [r for r in state.toxicity_risks if r.severity == "High"]
        tox_summary = f"고위험 {len(high_risks)}건: " + ", ".join([r.risk_type for r in high_risks]) if high_risks else "고위험 없음"
    
    patent_summary = "분석 실패"
    if state.patent_landscape:
        expired = [p for p in state.patent_landscape if p.status == "Expired"]
        active = [p for p in state.patent_landscape if p.status == "Active"]
        patent_summary = f"만료 특허 {len(expired)}건, 활성 특허 {len(active)}건"
    
    competitor_summary = "분석 실패"
    if state.competitors:
        approved = [c for c in state.competitors if c.phase == "Approved"]
        competitor_summary = f"승인 제품 {len(approved)}건, 총 경쟁사 {len(state.competitors)}개"
    
    prompt = ChatPromptTemplate.from_messages([
        ("system", REPORT_SYSTEM_PROMPT),
        ("user", """Generate final report for the following ADC analysis:

Configuration:
- Target: {target_name}
- Antibody: {antibody_type}  
- Payload: {payload_id}
- Linker: {linker_id}
- DAR: {dar}

Analysis Results:
- Structure: {structure_summary}
- Toxicity: {tox_summary}
- Patent: {patent_summary}
- Competition: {competitor_summary}

Generate a comprehensive assessment with grade and recommendation.""")
    ])
    
    chain = prompt | llm
    
    response = await chain.ainvoke({
        "target_name": state.input.target_name,
        "antibody_type": state.input.antibody_type,
        "payload_id": state.input.payload_id,
        "linker_id": state.input.linker_id,
        "dar": state.input.dar,
        "structure_summary": structure_summary,
        "tox_summary": tox_summary,
        "patent_summary": patent_summary,
        "competitor_summary": competitor_summary
    })
    
    # TODO: JSON 파싱 및 PDF 생성
    # 현재는 Mock 데이터 반환
    
    # 점수 계산 (각 분석 결과 기반)
    efficacy_score = 85 if state.structure_analysis and state.structure_analysis.stability_score > 80 else 70
    toxicity_score = 60 if state.toxicity_risks and any(r.severity == "High" for r in state.toxicity_risks) else 80
    properties_score = 75
    patent_score = 90 if state.patent_landscape and any(p.status == "Expired" for p in state.patent_landscape) else 70
    market_score = 70 if state.competitors and len([c for c in state.competitors if c.phase == "Approved"]) > 2 else 80
    
    # 종합 등급 계산
    avg_score = (efficacy_score + toxicity_score + properties_score + patent_score + market_score) / 5
    if avg_score >= 85:
        grade = "A"
        recommendation = "Go"
    elif avg_score >= 75:
        grade = "B+"
        recommendation = "Conditional Go"
    elif avg_score >= 65:
        grade = "B"
        recommendation = "Conditional Go"
    elif avg_score >= 55:
        grade = "C"
        recommendation = "Proceed with Caution"
    else:
        grade = "D"
        recommendation = "No Go"
    
    return {
        "grade": grade,
        "recommendation": recommendation,
        "summary": f"{state.input.target_name} 타겟 ADC 분석 결과, 구조적 안정성은 양호하나 {state.input.payload_id.upper()} 기반 독성 리스크가 확인되었습니다. 특허 상황은 유리하며, 경쟁 환경 분석 결과 차별화 포인트를 확보할 필요가 있습니다. 추가 전임상 연구 후 임상 진입을 권고합니다.",
        "key_strengths": [
            f"구조 안정성 {efficacy_score}/100",
            "주요 특허 만료로 FTO 확보 가능",
            f"{state.input.antibody_type} 항체 플랫폼 검증됨"
        ],
        "key_risks": [
            r.risk_type for r in (state.toxicity_risks[:2] if state.toxicity_risks else [])
        ] or ["독성 데이터 부족"],
        "next_steps": [
            "추가 in-vitro 독성 연구",
            "GLP 전임상 독성시험",
            "CMC 개발 착수"
        ],
        "scores": {
            "efficacy": efficacy_score,
            "toxicity": toxicity_score,
            "properties": properties_score,
            "patent": patent_score,
            "market": market_score
        },
        "report_url": None  # TODO: PDF 생성 후 S3 URL
    }
