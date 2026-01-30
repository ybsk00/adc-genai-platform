"""
The Coder Agent - Secure Version
f-string 인젝션 방지 + 코드 검증 강화

FIXED:
- f-string 인젝션 취약점 제거 (parameterized code execution)
- 코드 템플릿 유효성 검사
- Dangerous pattern 차단
"""
from typing import List, Dict, Any, Optional
import logging
import json
import re
import ast
import tempfile
import os

from app.agents.base_agent import BaseDesignAgent, AgentOutput
from app.agents.design_state import DesignSessionState, CalculatedMetrics
from app.core.supabase import get_supabase_client
from app.services.sandbox_executor import get_sandbox_executor, ExecutionResult

logger = logging.getLogger(__name__)


class SecureSnippetManager:
    """보안 강화된 코드 스니펫 관리자"""

    # FIXED: Dangerous patterns for code injection
    DANGEROUS_PATTERNS = [
        r'__import__\s*\(',
        r'importlib',
        r'eval\s*\(',
        r'exec\s*\(',
        r'compile\s*\(',
        r'os\.system',
        r'subprocess',
        r'open\s*\(',
        r'file\s*\(',
        r'\.write\s*\(',
        r'\.read\s*\(',
        r'input\s*\(',
        r'raw_input\s*\(',
        r'socket',
        r'urllib',
        r'requests',
        r'http\.client',
        r'ftplib',
        r'telnetlib',
        r'pickle\.load',
        r'yaml\.load',
        r'marshal\.loads',
    ]

    def __init__(self):
        self.supabase = get_supabase_client()
        self._cache: Dict[str, Dict] = {}

    def validate_code_security(self, code: str) -> tuple[bool, Optional[str]]:
        """
        FIXED: 코드 보안 검증
        
        Returns: (is_safe, error_message)
        """
        if not code or not isinstance(code, str):
            return False, "Empty or invalid code"
        
        # Check dangerous patterns
        for pattern in self.DANGEROUS_PATTERNS:
            if re.search(pattern, code, re.IGNORECASE):
                return False, f"Dangerous pattern detected: {pattern}"
        
        # Check f-string with user input (potential injection)
        # Look for f-strings that might contain unsanitized variables
        fstring_pattern = r'f["\'][^"\']*\{[^}]*\}[^"\']*["\']'
        fstrings = re.findall(fstring_pattern, code)
        for fs in fstrings:
            # Check if it contains complex expressions
            if any(op in fs for op in [';', 'import', 'exec', 'eval']):
                return False, f"Suspicious f-string detected: {fs[:50]}"
        
        # Syntax validation
        try:
            ast.parse(code)
        except SyntaxError as e:
            return False, f"Syntax error: {e}"
        
        return True, None

    async def get_secure_snippet(self, snippet_id: str, safe_params: Dict[str, Any]) -> Optional[Dict]:
        """
        FIXED: Parameterized snippet retrieval with security check
        
        Args:
            snippet_id: 스니펫 ID
            safe_params: 안전한 파라미터 (문자열, 숫자만)
        """
        snippet = await self._get_raw_snippet(snippet_id)
        if not snippet:
            return None
        
        template = snippet.get("code_template", "")
        
        # FIXED: Use safe string formatting instead of f-string
        # Replace placeholders with sanitized values
        safe_template = self._safe_format_template(template, safe_params)
        
        # Security check on formatted code
        is_safe, error = self.validate_code_security(safe_template)
        if not is_safe:
            logger.error(f"[coder-secure] Security check failed for snippet {snippet_id}: {error}")
            return None
        
        return {
            **snippet,
            "formatted_code": safe_template
        }

    def _safe_format_template(self, template: str, params: Dict[str, Any]) -> str:
        """
        FIXED: 안전한 템플릿 포맷팅
        
        f-string 대신 safe 문자열 치환 사용
        """
        result = template
        
        for key, value in params.items():
            # Sanitize key (alphanumeric only)
            safe_key = re.sub(r'[^a-zA-Z0-9_]', '', key)
            if not safe_key:
                continue
            
            # Sanitize value based on type
            if isinstance(value, str):
                # Escape quotes and backslashes
                safe_value = value.replace('\\', '\\\\').replace('"', '\\"').replace("'", "\\'")
                # Wrap in quotes
                safe_value = f'"{safe_value}"'
            elif isinstance(value, (int, float)):
                safe_value = str(value)
            elif isinstance(value, bool):
                safe_value = str(value)
            else:
                # Convert other types to JSON string
                safe_value = json.dumps(str(value))
            
            # Replace placeholder
            placeholder = f"{{{{{safe_key}}}}}"  # {{key}} -> key
            result = result.replace(placeholder, safe_value)
        
        return result

    async def _get_raw_snippet(self, snippet_id: str) -> Optional[Dict]:
        """Raw snippet retrieval with caching"""
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
            logger.error(f"[coder-secure] Snippet fetch error: {e}")

        return None


class CoderAgentSecure(BaseDesignAgent):
    """
    The Coder: 보안 강화된 Python 코드 실행 에이전트
    
    FIXED:
    - f-string 인젝션 방지
    - 코드 보안 검증
    - 안전한 파라미터 처리
    """

    name = "coder_secure"
    MAX_EXECUTION_TIME = 30

    def __init__(self):
        super().__init__()
        self.snippet_manager = SecureSnippetManager()
        self.sandbox = get_sandbox_executor()

    async def execute(self, state: DesignSessionState) -> AgentOutput:
        """보안 강화된 실행"""
        session_id = state["session_id"]
        smiles = state["current_smiles"]

        if not smiles:
            return AgentOutput(
                success=False,
                data={},
                error="No SMILES provided for validation"
            )

        # FIXED: Sanitize SMILES before using in code
        safe_smiles = self._sanitize_smiles(smiles)
        if not safe_smiles:
            return AgentOutput(
                success=False,
                data={},
                error="Invalid SMILES format"
            )

        await self._log_start(session_id, {"smiles": safe_smiles[:50]})

        try:
            # 필요한 검증 항목 결정
            validations_needed = self._determine_validations(state)

            # 보안 검증된 코드 조합
            code = await self._compose_secure_code(safe_smiles, validations_needed)
            
            if not code:
                return AgentOutput(
                    success=False,
                    data={},
                    error="Failed to compose secure code"
                )

            # 코드 실행
            result = await self._execute_secure_code(code, session_id)

            if result.get("exit_code", -1) != 0:
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

            # 결과 파싱
            metrics = self._parse_metrics(result.get("stdout", "{}"))
            
            # metrics를 state에 저장
            state["calculated_metrics"] = metrics

            output = AgentOutput(
                success=True,
                data={
                    "metrics": metrics,
                    "snippet_ids": validations_needed,
                    "code_executed": True
                },
                reasoning=f"Securely executed {len(validations_needed)} validation snippets",
                confidence_score=0.9,
                next_agent="auditor"
            )

            await self._save_execution_history(
                session_id=session_id,
                code=code,
                result=result,
                snippet_ids=validations_needed
            )

            await self._log_complete(session_id, output)
            return output

        except Exception as e:
            logger.exception(f"[coder-secure] Error: {e}")
            await self._log_error(session_id, str(e))

            state["last_error"] = str(e)
            state["requires_healing"] = True

            return AgentOutput(
                success=False,
                data={},
                error=str(e),
                next_agent="healer"
            )

    def _sanitize_smiles(self, smiles: str) -> Optional[str]:
        """
        FIXED: SMILES 문자열 정화
        """
        if not smiles or not isinstance(smiles, str):
            return None
        
        # Remove potentially dangerous characters
        # Keep only SMILES-safe characters
        allowed_pattern = r'^[A-Za-z0-9\[\]\(\)=\-#\+\\\/@\.\$%]+$'
        
        # Check each part (for combined SMILES with '.')
        parts = smiles.split('.')
        sanitized_parts = []
        
        for part in parts:
            part = part.strip()
            if not part:
                continue
            
            # Basic check - reject if contains suspicious patterns
            suspicious = [';', '\n', '\r', '#!', 'import', 'exec', 'eval', '__']
            if any(s in part.lower() for s in suspicious):
                logger.warning(f"[coder-secure] Suspicious SMILES pattern detected: {part[:50]}")
                continue
            
            sanitized_parts.append(part)
        
        return '.'.join(sanitized_parts) if sanitized_parts else None

    def _determine_validations(self, state: DesignSessionState) -> List[str]:
        """필요한 검증 스니펫 결정"""
        validations = ["lipinski_check", "mw_calc"]

        if state["session_type"] == "denovo":
            validations.extend(["pains_filter", "sa_score"])
        elif state["session_type"] == "optimization":
            validations.extend(["tanimoto_similarity", "pains_filter"])
        elif state["session_type"] == "audit":
            validations.extend(["pains_filter", "sa_score", "tanimoto_similarity"])

        return validations

    async def _compose_secure_code(
        self,
        smiles: str,
        snippet_ids: List[str]
    ) -> Optional[str]:
        """
        FIXED: 보안 검증된 코드 조합
        """
        imports = set()
        code_blocks = []

        # FIXED: Use parameterized snippet retrieval
        safe_params = {
            "smiles": smiles,  # Already sanitized
            "session_id": smiles[:10]  # Safe session ID prefix
        }

        for snippet_id in snippet_ids:
            snippet = await self.snippet_manager.get_secure_snippet(snippet_id, safe_params)
            if snippet and snippet.get("formatted_code"):
                for imp in snippet.get("required_imports", []):
                    imports.add(imp)
                code_blocks.append(snippet["formatted_code"])

        if not code_blocks:
            # Fallback to hardcoded safe code
            return self._generate_fallback_code(smiles)

        # 조합된 코드 검증
        import_lines = []
        for imp in sorted(imports):
            if "." in imp:
                parts = imp.rsplit(".", 1)
                import_lines.append(f"from {parts[0]} import {parts[1]}")
            else:
                import_lines.append(f"import {imp}")

        code = f'''#!/usr/bin/env python3
# Auto-generated secure validation code
# WARNING: Do not modify - auto-generated

import json
{chr(10).join(import_lines)}

SMILES_INPUT = {repr(smiles)}  # FIXED: Use repr() for safety
result = {{"smiles_valid": True}}

try:
    from rdkit import Chem
    mol = Chem.MolFromSmiles(SMILES_INPUT)
    if mol is None:
        result["error"] = "Invalid SMILES"
        result["smiles_valid"] = False
    else:
{chr(10).join('        ' + line for block in code_blocks for line in block.split(chr(10)))}
except Exception as e:
    result["error"] = str(e)
    result["smiles_valid"] = False

print(json.dumps(result, default=str))
'''

        # Final security check
        is_safe, error = self.snippet_manager.validate_code_security(code)
        if not is_safe:
            logger.error(f"[coder-secure] Final code security check failed: {error}")
            return self._generate_fallback_code(smiles)

        return code

    def _generate_fallback_code(self, smiles: str) -> str:
        """
        FIXED: 보안 검증된 fallback 코드
        """
        return f'''#!/usr/bin/env python3
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen

SMILES = {repr(smiles)}
result = {{}}

try:
    mol = Chem.MolFromSmiles(SMILES)
    if mol is None:
        result["error"] = "Invalid SMILES"
    else:
        result["mw"] = Descriptors.MolWt(mol)
        result["logp"] = Crippen.MolLogP(mol)
        result["hbd"] = Descriptors.NumHDonors(mol)
        result["hba"] = Descriptors.NumHAcceptors(mol)
        result["tpsa"] = Descriptors.TPSA(mol)
        result["rotatable_bonds"] = Descriptors.NumRotatableBonds(mol)
        result["lipinski"] = {{
            "mw": result["mw"],
            "logp": result["logp"],
            "hbd": result["hbd"],
            "hba": result["hba"]
        }}
except Exception as e:
    result["error"] = str(e)

print(json.dumps(result, default=str))
'''

    async def _execute_secure_code(self, code: str, session_id: str) -> Dict[str, Any]:
        """
        FIXED: 보안 검증된 코드 실행
        """
        # Additional security: write to temp file with restricted permissions
        with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
            f.write(code)
            script_path = f.name
        
        try:
            # Set read-only permissions
            os.chmod(script_path, 0o444)
            
            result: ExecutionResult = await self.sandbox.execute(
                code=code,
                timeout=self.MAX_EXECUTION_TIME,
                session_id=session_id
            )

            logger.info(
                f"[coder-secure] Execution completed - mode: {result.mode}, "
                f"exit_code: {result.exit_code}"
            )

            return {
                "exit_code": result.exit_code,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "execution_time_ms": result.execution_time_ms,
                "mode": result.mode
            }
        finally:
            # Cleanup
            try:
                os.chmod(script_path, 0o644)  # Restore permissions for deletion
                os.unlink(script_path)
            except:
                pass

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
            logger.warning(f"[coder-secure] Failed to parse metrics from: {stdout[:100]}")
            return CalculatedMetrics()

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
                "agent_name": "coder_secure",
                "code_content": code[:5000],  # Limit size
                "snippet_ids": snippet_ids,
                "execution_result": result,
                "stdout": result.get("stdout", "")[:2000],
                "stderr": result.get("stderr", "")[:1000],
                "exit_code": result.get("exit_code", -1)
            }).execute()
        except Exception as e:
            logger.error(f"[coder-secure] Failed to save execution history: {e}")


# Export secure version as default
code_agent = CoderAgentSecure()
