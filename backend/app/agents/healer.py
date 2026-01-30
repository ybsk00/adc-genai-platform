"""
The Healer Agent - 자가 치유
코드 에러 자동 수정 및 Knowledge Loop
"""
from typing import List, Dict, Any, Optional
import logging
import re
from datetime import datetime

from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.messages import HumanMessage, SystemMessage

from app.agents.base_agent import BaseDesignAgent, AgentOutput
from app.agents.design_state import DesignSessionState
from app.core.config import settings
from app.core.supabase import get_supabase_client

logger = logging.getLogger(__name__)


class HealerAgent(BaseDesignAgent):
    """
    The Healer: 자가 치유 에이전트

    핵심 기능:
    1. 에러 유형 분류
    2. LLM 기반 코드 자동 수정
    3. 최대 3회 재시도
    4. Knowledge Loop: 성공적 수정 패턴 자동 등록

    Constraint 2: 수정 횟수와 에러 메시지를 기록하여 모델 튜닝 데이터로 활용
    Constraint 5: 성공적 수정 패턴을 candidate_snippets에 자동 등록
    """

    name = "healer"
    MAX_RETRIES = 3

    # 에러 유형 분류
    ERROR_TYPES = {
        "SYNTAX_ERROR": ["syntaxerror", "indentationerror"],
        "CHEMICAL_VALIDITY_ERROR": ["valence", "kekulization", "sanitization", "atom"],
        "IMPORT_ERROR": ["importerror", "modulenotfounderror", "no module"],
        "TIMEOUT_ERROR": ["timeout", "timed out"],
        "MEMORY_ERROR": ["memoryerror", "killed"],
        "TYPE_ERROR": ["typeerror", "cannot convert", "expected"],
        "VALUE_ERROR": ["valueerror", "invalid"],
    }

    def __init__(self):
        super().__init__()
        self.llm = ChatGoogleGenerativeAI(
            model="gemini-2.0-flash-exp",
            google_api_key=settings.GOOGLE_API_KEY,
            temperature=0.2
        )

    async def execute(self, state: DesignSessionState) -> AgentOutput:
        """메인 실행"""
        session_id = state["session_id"]
        original_code = state.get("last_code", "")
        error_message = state.get("last_error", "")
        retry_count = state.get("healing_attempts", 0)

        # 재시도 횟수 초과
        if retry_count >= self.MAX_RETRIES:
            await self._log_healing(session_id, retry_count, [], False)
            return AgentOutput(
                success=False,
                data={
                    "status": "manual_review_required",
                    "total_attempts": retry_count,
                    "final_error": error_message
                },
                reasoning="Maximum retry attempts (3) exceeded. Manual review required.",
                error="Max healing retries exceeded"
            )

        await self._log_start(session_id, {
            "retry_count": retry_count + 1,
            "error_preview": error_message[:200] if error_message else None
        })

        try:
            # 1. 기존 수정 패턴 검색 (Knowledge Loop)
            existing_fix = await self._search_candidate_snippets(error_message)

            if existing_fix:
                logger.info(f"[healer] Found existing fix pattern: {existing_fix['id']}")
                fixed_code = self._apply_existing_fix(original_code, existing_fix)
                fix_explanation = f"Applied existing fix pattern: {existing_fix['fix_strategy']}"
            else:
                # 2. 에러 유형 분류
                error_type = self._classify_error(error_message)

                # 3. LLM을 사용하여 수정된 코드 생성
                fixed_code, fix_explanation = await self._generate_fix(
                    original_code=original_code,
                    error_message=error_message,
                    error_type=error_type,
                    attempt=retry_count + 1
                )

            # 4. 수정 내역 로깅 (모델 튜닝용)
            fix_log = {
                "attempt": retry_count + 1,
                "error_type": self._classify_error(error_message),
                "original_error": error_message[:500],
                "fix_explanation": fix_explanation,
                "used_existing_pattern": existing_fix is not None,
                "timestamp": datetime.utcnow().isoformat()
            }

            await self._save_fix_log(session_id, fix_log)

            output = AgentOutput(
                success=True,
                data={
                    "fixed_code": fixed_code,
                    "retry_count": retry_count + 1,
                    "fix_log": fix_log
                },
                reasoning=f"Error type: {fix_log['error_type']}. Applied fix: {fix_explanation}",
                confidence_score=0.7 if existing_fix else 0.5,
                next_agent="coder"
            )

            await self._log_complete(session_id, output)
            return output

        except Exception as e:
            logger.exception(f"[healer] Error: {e}")
            await self._log_error(session_id, str(e))
            return AgentOutput(
                success=False,
                data={},
                error=str(e)
            )

    def _classify_error(self, error_message: str) -> str:
        """에러 유형 분류"""
        error_lower = error_message.lower()

        for error_type, keywords in self.ERROR_TYPES.items():
            for keyword in keywords:
                if keyword in error_lower:
                    return error_type

        return "RUNTIME_ERROR"

    async def _search_candidate_snippets(self, error_message: str) -> Optional[Dict]:
        """
        Knowledge Loop: 기존 성공 패턴 검색

        candidate_snippets 테이블에서 유사 에러 패턴 검색
        """
        error_type = self._classify_error(error_message)

        try:
            result = self.supabase.table("candidate_snippets").select("*").eq(
                "error_type", error_type
            ).eq(
                "status", "candidate"  # 또는 'promoted'
            ).gte(
                "confidence_score", 0.7
            ).order(
                "success_count", desc=True
            ).limit(1).execute()

            if result.data:
                return result.data[0]

        except Exception as e:
            logger.warning(f"[healer] Candidate snippet search error: {e}")

        return None

    def _apply_existing_fix(self, original_code: str, fix_pattern: Dict) -> str:
        """기존 수정 패턴 적용"""
        fix_code = fix_pattern.get("fix_code", "")
        fix_strategy = fix_pattern.get("fix_strategy", "")

        # 전략에 따라 수정 적용
        if fix_strategy == "add_null_check":
            # mol is None 체크 추가
            if "mol = Chem.MolFromSmiles" in original_code and "if mol is None" not in original_code:
                original_code = original_code.replace(
                    "mol = Chem.MolFromSmiles(smiles)",
                    "mol = Chem.MolFromSmiles(smiles)\n    if mol is None:\n        result['error'] = 'Invalid SMILES'\n    else:"
                )
        elif fix_strategy == "add_exception_handler":
            # try-except 추가
            if "try:" not in original_code:
                original_code = f"try:\n{original_code}\nexcept Exception as e:\n    result['error'] = str(e)"

        return original_code

    async def _generate_fix(
        self,
        original_code: str,
        error_message: str,
        error_type: str,
        attempt: int
    ) -> tuple[str, str]:
        """LLM을 사용하여 수정된 코드 생성"""

        system_prompt = """You are an expert Python developer specializing in cheminformatics (RDKit).
Your task is to fix Python code errors, especially those related to molecular processing.

RULES:
1. Only modify what's necessary to fix the error
2. Keep the code structure intact
3. Add proper error handling
4. Ensure RDKit functions handle None/invalid molecules
5. Return valid Python code

Response format:
FIXED_CODE:
```python
<your fixed code here>
```

EXPLANATION:
<brief explanation of what you fixed>"""

        user_prompt = f"""Fix this Python code error:

ERROR TYPE: {error_type}
ERROR MESSAGE:
{error_message}

ORIGINAL CODE:
```python
{original_code}
```

ATTEMPT: {attempt} of 3

Please fix the code and explain what you changed."""

        try:
            response = await self.llm.ainvoke([
                SystemMessage(content=system_prompt),
                HumanMessage(content=user_prompt)
            ])

            content = response.content

            # Parse fixed code
            fixed_code = original_code
            explanation = "Unable to extract fix"

            if "```python" in content:
                code_match = re.search(r'```python\n(.*?)```', content, re.DOTALL)
                if code_match:
                    fixed_code = code_match.group(1).strip()

            if "EXPLANATION:" in content:
                explanation = content.split("EXPLANATION:")[-1].strip()
            elif "fixed" in content.lower():
                # 간단한 설명 추출
                explanation = content.split("\n")[-1][:200]

            return fixed_code, explanation

        except Exception as e:
            logger.error(f"[healer] LLM fix generation error: {e}")
            # Fallback: 기본 에러 핸들링 추가
            fallback_code = self._add_basic_error_handling(original_code)
            return fallback_code, f"Added basic error handling (LLM failed: {str(e)[:50]})"

    def _add_basic_error_handling(self, code: str) -> str:
        """기본 에러 핸들링 추가"""
        if "try:" not in code:
            lines = code.split("\n")
            indented_lines = ["    " + line for line in lines]
            return f"""try:
{chr(10).join(indented_lines)}
except Exception as e:
    result["error"] = str(e)
    print(json.dumps(result))"""
        return code

    async def _save_fix_log(self, session_id: str, fix_log: Dict):
        """
        Constraint 2: 수정 이력을 DB에 저장하여 모델 튜닝 데이터로 활용
        """
        try:
            self.supabase.table("agent_execution_logs").insert({
                "session_id": session_id,
                "agent_name": "healer",
                "status": "healed",
                "reasoning": fix_log["fix_explanation"],
                "fix_logs": [fix_log],
                "retry_count": fix_log["attempt"]
            }).execute()
        except Exception as e:
            logger.error(f"[healer] Failed to save fix log: {e}")

    async def register_successful_fix(
        self,
        session_id: str,
        error_message: str,
        fixed_code: str,
        fix_explanation: str
    ):
        """
        Knowledge Loop (Constraint 5): 성공적인 수정 패턴을 candidate_snippets에 등록

        Coder가 수정된 코드 실행 성공 시 호출
        """
        error_type = self._classify_error(error_message)
        error_pattern = self._extract_error_pattern(error_message)
        fix_strategy = self._classify_fix_strategy(fix_explanation)

        try:
            # 기존 패턴 확인
            existing = self.supabase.table("candidate_snippets").select("id, success_count").eq(
                "error_type", error_type
            ).eq(
                "fix_strategy", fix_strategy
            ).single().execute()

            if existing.data:
                # 기존 패턴 성공 카운트 증가
                self.supabase.table("candidate_snippets").update({
                    "success_count": existing.data["success_count"] + 1,
                    "confidence_score": min(0.95, 0.5 + (existing.data["success_count"] + 1) * 0.05),
                    "updated_at": datetime.utcnow().isoformat()
                }).eq("id", existing.data["id"]).execute()
            else:
                # 새 패턴 등록
                self.supabase.table("candidate_snippets").insert({
                    "source_session_id": session_id,
                    "error_type": error_type,
                    "error_pattern": error_pattern,
                    "error_context": {"full_error": error_message[:1000]},
                    "fix_code": fixed_code,
                    "fix_explanation": fix_explanation,
                    "fix_strategy": fix_strategy,
                    "confidence_score": 0.50
                }).execute()

            logger.info(f"[healer] Registered successful fix pattern: {fix_strategy}")

        except Exception as e:
            logger.error(f"[healer] Failed to register fix pattern: {e}")

    def _extract_error_pattern(self, error_message: str) -> str:
        """에러 메시지에서 재사용 가능한 패턴 추출"""
        # 파일 경로 및 라인 번호 제거
        pattern = re.sub(r'File ".*?", line \d+', 'File "...", line X', error_message)
        # 변수명 제거
        pattern = re.sub(r"'[^']*'", "'VAR'", pattern)
        return pattern[:500]

    def _classify_fix_strategy(self, fix_explanation: str) -> str:
        """수정 전략 분류"""
        explanation = fix_explanation.lower()

        if "null check" in explanation or "none check" in explanation:
            return "add_null_check"
        elif "try" in explanation or "except" in explanation:
            return "add_exception_handler"
        elif "import" in explanation:
            return "fix_import"
        elif "type" in explanation or "cast" in explanation:
            return "fix_type_conversion"
        elif "syntax" in explanation:
            return "fix_syntax"
        elif "indent" in explanation:
            return "fix_indentation"
        else:
            return "general_fix"
