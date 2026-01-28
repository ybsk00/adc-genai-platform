"""
The Alchemist Agent - ADC 후보 물질 설계
Golden Set 기반 SMILES 생성 및 최적화
"""
from typing import List, Dict, Any, Optional
import logging
import json

from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.messages import HumanMessage, SystemMessage

from app.agents.base_agent import BaseDesignAgent, AgentOutput
from app.agents.design_state import DesignSessionState, CandidateStructure
from app.core.config import settings
from app.core.supabase import get_supabase_client

logger = logging.getLogger(__name__)


class AlchemistAgent(BaseDesignAgent):
    """
    The Alchemist: SMILES 설계 에이전트

    핵심 기능:
    1. Golden Set에서 유사 구조 검색
    2. 타겟/적응증에 맞는 후보 SMILES 생성
    3. DAR, Linker 선호도 반영
    4. 상위 N개 후보 반환
    """

    name = "alchemist"

    def __init__(self):
        super().__init__()
        self.llm = ChatGoogleGenerativeAI(
            model="gemini-2.0-flash-exp",
            google_api_key=settings.GEMINI_API_KEY,
            temperature=0.3
        )

    async def execute(self, state: DesignSessionState) -> AgentOutput:
        """메인 실행"""
        session_id = state["session_id"]
        await self._log_start(session_id, {
            "target": state["target_antigen"],
            "dar": state["requested_dar"],
            "linker": state["linker_preference"]
        })

        try:
            # 1. Golden Set에서 관련 구조 검색
            golden_set_refs = await self._search_golden_set(
                target=state["target_antigen"],
                indication=state["target_indication"],
                limit=10
            )

            # 2. LLM을 사용하여 후보 SMILES 생성
            candidates = await self._generate_candidates(
                state=state,
                golden_set_refs=golden_set_refs
            )

            # 3. 결과 구성
            primary_smiles = candidates[0]["smiles"] if candidates else ""

            output = AgentOutput(
                success=True,
                data={
                    "primary_smiles": primary_smiles,
                    "candidates": candidates,
                    "golden_set_refs": [g["id"] for g in golden_set_refs],
                    "summary": f"Generated {len(candidates)} candidates based on {len(golden_set_refs)} Golden Set references"
                },
                reasoning=f"Target: {state['target_antigen']}, DAR: {state['requested_dar']}, "
                         f"Found {len(golden_set_refs)} similar structures in Golden Set",
                confidence_score=0.85 if golden_set_refs else 0.5,
                next_agent="coder"
            )

            await self._log_complete(session_id, output)
            return output

        except Exception as e:
            logger.exception(f"[alchemist] Error: {e}")
            await self._log_error(session_id, str(e))
            return AgentOutput(
                success=False,
                data={},
                error=str(e)
            )

    async def _search_golden_set(
        self,
        target: str,
        indication: str,
        limit: int = 10
    ) -> List[Dict[str, Any]]:
        """
        Golden Set에서 유사 구조 검색

        타겟 항원과 적응증 기반 검색
        """
        try:
            # 타겟 기반 검색
            query = self.supabase.table("golden_set").select(
                "id, drug_name, target_1, target_2, indication, "
                "payload_smiles, linker_smiles, linker_type, dar, "
                "clinical_status, orr_pct, os_months"
            )

            if target:
                query = query.or_(f"target_1.ilike.%{target}%,target_2.ilike.%{target}%")

            result = query.limit(limit).execute()

            if result.data:
                logger.info(f"[alchemist] Found {len(result.data)} Golden Set entries for target: {target}")
                return result.data

            # 타겟으로 못 찾으면 적응증으로 검색
            if indication:
                result = self.supabase.table("golden_set").select("*").ilike(
                    "indication", f"%{indication}%"
                ).limit(limit).execute()
                return result.data or []

            return []

        except Exception as e:
            logger.error(f"[alchemist] Golden Set search error: {e}")
            return []

    async def _generate_candidates(
        self,
        state: DesignSessionState,
        golden_set_refs: List[Dict]
    ) -> List[CandidateStructure]:
        """
        LLM을 사용하여 후보 SMILES 생성

        Golden Set 참조와 사용자 요구사항을 기반으로 생성
        """
        # Golden Set 정보 요약
        ref_summary = self._summarize_golden_set(golden_set_refs)

        system_prompt = """You are an expert medicinal chemist specializing in Antibody-Drug Conjugates (ADCs).
Your task is to design ADC payload-linker combinations based on FDA-approved ADC structures (Golden Set).

IMPORTANT RULES:
1. Generate valid SMILES strings that can be parsed by RDKit
2. Consider the target antigen and indication
3. Respect the requested DAR (Drug-Antibody Ratio)
4. Consider linker preference (cleavable vs non-cleavable)
5. Prioritize structures similar to successful FDA-approved ADCs

Return your response in JSON format:
{
    "candidates": [
        {
            "smiles": "VALID_SMILES_STRING",
            "name": "Descriptive name",
            "rationale": "Why this structure was chosen",
            "similarity_to_ref": "Which Golden Set drug it's based on",
            "estimated_score": 0.0 to 1.0
        }
    ]
}"""

        user_prompt = f"""Design ADC payload-linker candidates with these requirements:

TARGET ANTIGEN: {state['target_antigen']}
INDICATION: {state['target_indication']}
DESIRED DAR: {state['requested_dar']}
LINKER PREFERENCE: {state['linker_preference']}
DESIGN GOAL: {state['design_goal']}

REFERENCE STRUCTURES (FDA-approved ADCs):
{ref_summary}

Generate 3-5 candidate structures ranked by predicted efficacy.
Focus on structures similar to successful references but optimized for the target."""

        try:
            response = await self.llm.ainvoke([
                SystemMessage(content=system_prompt),
                HumanMessage(content=user_prompt)
            ])

            # Parse JSON response
            content = response.content
            # Extract JSON from response
            if "```json" in content:
                content = content.split("```json")[1].split("```")[0]
            elif "```" in content:
                content = content.split("```")[1].split("```")[0]

            data = json.loads(content)

            candidates = []
            for idx, c in enumerate(data.get("candidates", [])):
                candidates.append(CandidateStructure(
                    smiles=c.get("smiles", ""),
                    rank=idx + 1,
                    score=float(c.get("estimated_score", 0.5)),
                    metrics=None,  # Coder가 계산
                    validation=None,  # Auditor가 검증
                    is_masked=state["tier"] == "free" and idx > 0
                ))

            return candidates

        except json.JSONDecodeError as e:
            logger.error(f"[alchemist] JSON parse error: {e}")
            # Fallback: Golden Set에서 직접 가져오기
            return self._fallback_candidates(golden_set_refs, state["tier"])

        except Exception as e:
            logger.error(f"[alchemist] LLM error: {e}")
            return self._fallback_candidates(golden_set_refs, state["tier"])

    def _summarize_golden_set(self, refs: List[Dict]) -> str:
        """Golden Set 참조 요약"""
        if not refs:
            return "No reference structures available."

        lines = []
        for ref in refs[:5]:  # 최대 5개
            lines.append(
                f"- {ref.get('drug_name', 'Unknown')}: "
                f"Target={ref.get('target_1', 'N/A')}, "
                f"DAR={ref.get('dar', 'N/A')}, "
                f"Linker={ref.get('linker_type', 'N/A')}, "
                f"ORR={ref.get('orr_pct', 'N/A')}%"
            )
        return "\n".join(lines)

    def _fallback_candidates(
        self,
        golden_set_refs: List[Dict],
        tier: str
    ) -> List[CandidateStructure]:
        """Fallback: Golden Set에서 직접 후보 생성"""
        candidates = []
        for idx, ref in enumerate(golden_set_refs[:3]):
            smiles = ref.get("payload_smiles") or ref.get("linker_smiles") or ""
            if smiles:
                candidates.append(CandidateStructure(
                    smiles=smiles,
                    rank=idx + 1,
                    score=0.7,
                    metrics=None,
                    validation=None,
                    is_masked=tier == "free" and idx > 0
                ))
        return candidates
