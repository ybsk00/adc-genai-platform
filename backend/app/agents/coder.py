"""
The Coder Agent - Python 코드 실행
검증된 Snippet Library를 조합하여 Sandbox에서 실행

Constraint 7: Sandbox Library Version Lock
"""
from typing import List, Dict, Any, Optional
import logging
import json

from app.agents.base_agent import BaseDesignAgent, AgentOutput
from app.agents.design_state import DesignSessionState, CalculatedMetrics
from app.core.supabase import get_supabase_client
from app.services.sandbox_executor import get_sandbox_executor, ExecutionResult

logger = logging.getLogger(__name__)


class SnippetManager:
    """
    검증된 코드 스니펫 관리자

    DB에서 스니펫을 로드하고 템플릿을 채워 실행 가능한 코드 생성
    """

    def __init__(self):
        self.supabase = get_supabase_client()
        self._cache: Dict[str, Dict] = {}

    async def get(self, snippet_id: str) -> Optional[Dict]:
        """스니펫 조회 (캐시 사용)"""
        if snippet_id in self._cache:
            return self._cache[snippet_id]

        try:
            result = self.supabase.table("code_snippet_library").select("*").eq(
                "id", snippet_id
            ).single().execute()

            if result.data:
                self._cache[snippet_id] = result.data
                return result.data
        except Exception as e:
            logger.error(f"[coder] Snippet fetch error: {e}")

        return None

    async def get_all_verified(self) -> List[Dict]:
        """검증된 스니펫 전체 조회"""
        try:
            result = self.supabase.table("code_snippet_library").select("*").eq(
                "is_verified", True
            ).execute()
            return result.data or []
        except Exception as e:
            logger.error(f"[coder] Verified snippets fetch error: {e}")
            return []


class CoderAgent(BaseDesignAgent):
    """
    The Coder: Python 코드 실행 에이전트

    핵심 기능:
    1. 검증된 Snippet Library에서 코드 조합
    2. Docker Sandbox에서 안전하게 실행
    3. 에러 발생 시 The Healer로 전달
    4. 실행 결과로 물성 계산
    """

    name = "coder"
    MAX_EXECUTION_TIME = 30  # seconds

    def __init__(self):
        super().__init__()
        self.snippet_manager = SnippetManager()
        self.sandbox = get_sandbox_executor()

    async def execute(self, state: DesignSessionState) -> AgentOutput:
        """메인 실행"""
        session_id = state["session_id"]
        smiles = state["current_smiles"]

        if not smiles:
            return AgentOutput(
                success=False,
                data={},
                error="No SMILES provided for validation"
            )

        await self._log_start(session_id, {"smiles": smiles[:50] + "..." if len(smiles) > 50 else smiles})

        try:
            # 1. 필요한 검증 항목 결정
            validations_needed = self._determine_validations(state)
            logger.info(f"[coder] session_type={state['session_type']}, validations={validations_needed}")

            # 2. 검증된 스니펫으로 코드 조합
            code = await self._compose_code_from_snippets(smiles, validations_needed)
            logger.info(f"[coder] Generated code length={len(code)} chars, first 200: {code[:200]}")

            # 3. 코드 실행 (현재는 직접 실행, 추후 Docker Sandbox로 대체)
            result = await self._execute_code(code, session_id)
            logger.info(f"[coder] Execution result: exit_code={result.get('exit_code')}, "
                        f"stdout_len={len(result.get('stdout', ''))}, "
                        f"stderr={result.get('stderr', '')[:300]}")

            if result.get("exit_code", -1) != 0:
                # Error -> The Healer로 전달
                state["last_code"] = code
                state["last_error"] = result.get("stderr", "Unknown error")
                state["requires_healing"] = True

                await self._log_error(session_id, result.get("stderr", ""))

                return AgentOutput(
                    success=False,
                    data={
                        "code": code,
                        "error": result.get("stderr"),
                        "snippet_ids": validations_needed
                    },
                    error=result.get("stderr"),
                    next_agent="healer"
                )

            # 4. 결과 파싱 및 메트릭 추출
            metrics = self._parse_metrics(result.get("stdout", "{}"))

            output = AgentOutput(
                success=True,
                data={
                    "metrics": metrics,
                    "snippet_ids": validations_needed,
                    "code_executed": True
                },
                reasoning=f"Executed {len(validations_needed)} validation snippets successfully",
                confidence_score=0.9,
                next_agent="auditor"
            )

            # 코드 실행 이력 저장
            await self._save_execution_history(
                session_id=session_id,
                code=code,
                result=result,
                snippet_ids=validations_needed
            )

            await self._log_complete(session_id, output)
            return output

        except Exception as e:
            logger.exception(f"[coder] Error: {e}")
            await self._log_error(session_id, str(e))

            state["last_error"] = str(e)
            state["requires_healing"] = True

            return AgentOutput(
                success=False,
                data={},
                error=str(e),
                next_agent="healer"
            )

    def _determine_validations(self, state: DesignSessionState) -> List[str]:
        """필요한 검증 스니펫 결정"""
        validations = ["lipinski_check", "mw_calc"]

        # 세션 타입에 따라 추가 검증
        session_type = state["session_type"]
        if session_type in ("denovo", "navigator"):
            validations.extend(["pains_filter", "sa_score"])
        elif session_type == "optimization":
            validations.extend(["tanimoto_similarity", "pains_filter"])
        elif session_type == "audit":
            validations.extend(["pains_filter", "sa_score", "tanimoto_similarity"])

        return validations

    async def _compose_code_from_snippets(
        self,
        smiles: str,
        snippet_ids: List[str]
    ) -> str:
        """검증된 스니펫들을 조합하여 실행 가능한 코드 생성"""
        imports = set()
        code_blocks = []

        for snippet_id in snippet_ids:
            snippet = await self.snippet_manager.get(snippet_id)
            if snippet:
                # 필요한 import 수집
                for imp in snippet.get("required_imports", []):
                    imports.add(imp)

                # 코드 템플릿에서 placeholder 교체
                template = snippet.get("code_template", "")
                code_blocks.append(template)

        # 스니펫이 하나도 없으면 인라인 폴백 코드 사용
        if not code_blocks:
            logger.warning("[coder] No snippets found in DB — using inline fallback code")
            return self._get_inline_fallback_code(smiles)

        # 최종 코드 조합
        import_lines = []
        for imp in imports:
            if "." in imp:
                parts = imp.rsplit(".", 1)
                import_lines.append(f"from {parts[0]} import {parts[1]}")
            else:
                import_lines.append(f"import {imp}")

        # SMILES 안전하게 이스케이프 (backslash, quote 등)
        safe_smiles = smiles.replace("\\", "\\\\").replace('"', '\\"')

        # 코드 블록 조합 (들여쓰기 8칸 = else: 블록 내부)
        indented_blocks = chr(10).join(
            '        ' + line
            for block in code_blocks
            for line in block.split(chr(10))
        )

        code = (
            "import json\n"
            + chr(10).join(import_lines) + "\n"
            + "\n"
            + f'smiles = "{safe_smiles}"\n'
            + "result = {}\n"
            + "\n"
            + "try:\n"
            + "    from rdkit import Chem\n"
            + "    mol = Chem.MolFromSmiles(smiles)\n"
            + "    if mol is None:\n"
            + '        result["error"] = "Invalid SMILES"\n'
            + "    else:\n"
            + indented_blocks + "\n"
            + "except Exception as e:\n"
            + '    result["error"] = str(e)\n'
            + "\n"
            + "print(json.dumps(result))\n"
        )
        return code

    @staticmethod
    def _get_inline_fallback_code(smiles: str) -> str:
        """
        DB 스니펫 없을 때 인라인 폴백 코드

        실제 RDKit 연산만 수행 (추정값/가짜 데이터 없음 — Fail-Fast 원칙 준수)
        """
        safe_smiles = smiles.replace("\\", "\\\\").replace('"', '\\"')
        return f'''
import json

smiles = "{safe_smiles}"
result = {{}}

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Lipinski

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        result["error"] = "Invalid SMILES — RDKit could not parse"
    else:
        # Lipinski descriptors (실제 계산)
        result["lipinski"] = {{
            "mw": round(Descriptors.MolWt(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
            "hbd": Lipinski.NumHDonors(mol),
            "hba": Lipinski.NumHAcceptors(mol),
        }}
        result["mw"] = result["lipinski"]["mw"]
        result["tpsa"] = round(Descriptors.TPSA(mol), 2)
        result["rotatable_bonds"] = Lipinski.NumRotatableBonds(mol)

        # PAINS filter
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        catalog = FilterCatalog(params)
        pains_match = catalog.HasMatch(mol)
        result["pains_alert"] = pains_match

        # SA Score
        try:
            from rdkit.Contrib.SA_Score import sascorer
            result["sa_score"] = round(sascorer.calculateScore(mol), 2)
        except Exception:
            pass  # SA Score 모듈 없으면 skip (에러는 아님)

        result["source"] = "inline_fallback_rdkit"

except ImportError as e:
    result["error"] = f"RDKit not available: {{e}}"
except Exception as e:
    result["error"] = f"Calculation failed: {{e}}"

print(json.dumps(result))
'''

    async def _execute_code(self, code: str, session_id: str) -> Dict[str, Any]:
        """
        Sandbox에서 코드 실행

        - 개발 모드: subprocess
        - 프로덕션 모드: Docker 컨테이너 (Air-gapped)
        """
        result: ExecutionResult = await self.sandbox.execute(
            code=code,
            timeout=self.MAX_EXECUTION_TIME,
            session_id=session_id
        )

        logger.info(
            f"[coder] Execution completed - mode: {result.mode}, "
            f"exit_code: {result.exit_code}, time: {result.execution_time_ms}ms"
        )

        return {
            "exit_code": result.exit_code,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "execution_time_ms": result.execution_time_ms,
            "mode": result.mode
        }

    def _parse_metrics(self, stdout: str) -> CalculatedMetrics:
        """stdout에서 메트릭 파싱"""
        try:
            data = json.loads(stdout)

            lipinski = data.get("lipinski", {})

            return CalculatedMetrics(
                mw=lipinski.get("mw") or data.get("mw"),
                logp=lipinski.get("logp"),
                hbd=lipinski.get("hbd"),
                hba=lipinski.get("hba"),
                psa=data.get("psa"),
                tpsa=data.get("tpsa"),
                rotatable_bonds=data.get("rotatable_bonds"),
                sa_score=data.get("sa_score")
            )
        except json.JSONDecodeError:
            logger.warning(f"[coder] Failed to parse metrics from: {stdout[:100]}")
            return {}

    async def _save_execution_history(
        self,
        session_id: str,
        code: str,
        result: Dict,
        snippet_ids: List[str]
    ):
        """코드 실행 이력 저장"""
        try:
            self.supabase.table("code_execution_history").insert({
                "session_id": session_id,
                "agent_name": "coder",
                "code_content": code,
                "snippet_ids": snippet_ids,
                "execution_result": result,
                "stdout": result.get("stdout", ""),
                "stderr": result.get("stderr", ""),
                "exit_code": result.get("exit_code", -1)
            }).execute()
        except Exception as e:
            logger.error(f"[coder] Failed to save execution history: {e}")
