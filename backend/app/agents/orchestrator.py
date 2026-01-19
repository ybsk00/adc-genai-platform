"""
LangGraph Orchestrator - ADC 분석 워크플로우
6개 에이전트를 병렬/순차 실행하여 ADC 분석 수행
"""
from typing import Dict, Any, Literal
from datetime import datetime
import asyncio

from langgraph.graph import StateGraph, END
from langchain_openai import ChatOpenAI

from app.agents.state import ADCState, AgentStatus, JobStatus
from app.agents.structure_agent import run_structure_agent
from app.agents.tox_agent import run_tox_agent
from app.agents.patent_agent import run_patent_agent
from app.agents.competitor_agent import run_competitor_agent
from app.agents.clinical_agent import run_clinical_agent
from app.agents.report_agent import run_report_agent
from app.agents.commercial_agent import run_commercial_agent
from app.core.config import settings



def create_orchestrator() -> StateGraph:
    """
    LangGraph 오케스트레이터 생성
    
    워크플로우:
    1. Structure Agent (먼저 실행 - 3D 구조 분석)
    2. [병렬 실행] Toxicology, Patent, Competitor Agents
    3. Clinical Agent (앞선 결과 종합 후)
    4. Report Agent (최종 리포트 생성)
    """
    
    # StateGraph 생성
    workflow = StateGraph(ADCState)
    
    # 노드 추가
    workflow.add_node("structure", structure_node)
    workflow.add_node("parallel_analysis", parallel_analysis_node)
    workflow.add_node("clinical", clinical_node)
    workflow.add_node("report", report_node)
    
    # 엣지 정의 (실행 순서)
    workflow.set_entry_point("structure")
    workflow.add_edge("structure", "parallel_analysis")
    workflow.add_edge("parallel_analysis", "clinical")
    workflow.add_edge("clinical", "report")
    workflow.add_edge("report", END)
    
    return workflow.compile()


async def structure_node(state: ADCState) -> ADCState:
    """구조 분석 노드"""
    state.update_agent_status("structure", AgentStatus.RUNNING)
    
    try:
        result = await run_structure_agent(state)
        state.structure_analysis = result
        state.update_agent_status("structure", AgentStatus.DONE)
    except Exception as e:
        state.update_agent_status("structure", AgentStatus.ERROR, str(e))
    
    return state


async def parallel_analysis_node(state: ADCState) -> ADCState:
    """
    병렬 분석 노드
    Toxicology, Patent, Competitor 에이전트를 동시에 실행
    """
    # 각 에이전트를 Running으로 설정
    state.update_agent_status("toxicology", AgentStatus.RUNNING)
    state.update_agent_status("patent", AgentStatus.RUNNING)
    state.update_agent_status("competitor", AgentStatus.RUNNING)
    state.update_agent_status("commercial", AgentStatus.RUNNING)

    
    # 병렬 실행
    async def run_with_error_handling(agent_func, agent_id: str):
        try:
            return await agent_func(state)
        except Exception as e:
            state.update_agent_status(agent_id, AgentStatus.ERROR, str(e))
            return None
    
    results = await asyncio.gather(
        run_with_error_handling(run_tox_agent, "toxicology"),
        run_with_error_handling(run_patent_agent, "patent"),
        run_with_error_handling(run_competitor_agent, "competitor"),
        run_with_error_handling(run_commercial_agent, "commercial"),
        return_exceptions=True

    )
    
    # 결과 적용
    if results[0] and not isinstance(results[0], Exception):
        state.toxicity_risks = results[0]
        state.update_agent_status("toxicology", AgentStatus.DONE)
    
    if results[1] and not isinstance(results[1], Exception):
        state.patent_landscape = results[1]
        state.update_agent_status("patent", AgentStatus.DONE)
    
    if results[2] and not isinstance(results[2], Exception):
        state.competitors = results[2]
        state.update_agent_status("competitor", AgentStatus.DONE)

    if results[3] and not isinstance(results[3], Exception):
        state.commercial_feasibility = results[3]
        state.update_agent_status("commercial", AgentStatus.DONE)

    
    return state


async def clinical_node(state: ADCState) -> ADCState:
    """임상 기획 노드"""
    state.update_agent_status("clinical", AgentStatus.RUNNING)
    
    try:
        result = await run_clinical_agent(state)
        state.clinical_protocol = result
        state.update_agent_status("clinical", AgentStatus.DONE)
    except Exception as e:
        state.update_agent_status("clinical", AgentStatus.ERROR, str(e))
    
    return state


async def report_node(state: ADCState) -> ADCState:
    """리포트 생성 노드 (최종)"""
    state.update_agent_status("report", AgentStatus.RUNNING)
    
    try:
        result = await run_report_agent(state)
        
        # 최종 결과 설정
        state.final_grade = result.get("grade", "B")
        state.recommendation = result.get("recommendation", "Conditional Go")
        state.executive_summary = result.get("summary", "")
        state.report_url = result.get("report_url")
        
        # 점수 업데이트
        state.scores = result.get("scores", state.scores)
        
        state.update_agent_status("report", AgentStatus.DONE)
        
        # 전체 상태 판단
        error_agents = [a for a in state.agents.values() if a.status == AgentStatus.ERROR]
        if error_agents:
            state.status = JobStatus.PARTIAL
        else:
            state.status = JobStatus.COMPLETED
            
    except Exception as e:
        state.update_agent_status("report", AgentStatus.ERROR, str(e))
        state.status = JobStatus.FAILED
    
    return state


async def run_adc_analysis(input_data: dict, job_id: str, user_id: str) -> ADCState:
    """
    ADC 분석 실행 엔트리포인트
    
    Args:
        input_data: ADCInput 데이터
        job_id: 작업 ID
        user_id: 사용자 ID
    
    Returns:
        ADCState: 최종 분석 결과
    """
    from app.agents.state import ADCInput
    
    # 초기 상태 생성
    initial_state = ADCState(
        job_id=job_id,
        user_id=user_id,
        input=ADCInput(**input_data),
        status=JobStatus.RUNNING
    )
    
    # 오케스트레이터 실행
    orchestrator = create_orchestrator()
    final_state = await orchestrator.ainvoke(initial_state)
    
    return final_state
