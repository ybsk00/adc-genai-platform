"""
Sandbox Executor Service
Python 코드 안전 실행 서비스

Constraint 7: Sandbox Library Version Lock
- 개발 모드: subprocess (빠른 개발)
- 프로덕션 모드: Docker 컨테이너 (Air-gapped, 보안 강화)
"""
import os
import sys
import json
import asyncio
import tempfile
import subprocess
import logging
import platform
from typing import Dict, Any, Optional
from enum import Enum
from dataclasses import dataclass

logger = logging.getLogger(__name__)


class ExecutionMode(Enum):
    """실행 모드"""
    SUBPROCESS = "subprocess"  # 개발용: 직접 실행
    DOCKER = "docker"          # 프로덕션: Docker 격리


@dataclass
class ExecutionResult:
    """실행 결과"""
    success: bool
    exit_code: int
    stdout: str
    stderr: str
    execution_time_ms: int
    mode: str


class SandboxExecutor:
    """
    Python 코드 Sandbox 실행기

    보안 제약:
    - 네트워크 접근 차단 (Docker)
    - 메모리 제한 512MB (Docker)
    - 실행 시간 제한 30초
    - 허용된 라이브러리만 사용 가능
    """

    # 허용된 라이브러리 목록
    ALLOWED_IMPORTS = {
        "rdkit", "rdkit.Chem", "rdkit.Chem.Descriptors",
        "rdkit.Chem.Lipinski", "rdkit.Chem.AllChem",
        "rdkit.Chem.FilterCatalog", "rdkit.DataStructs",
        "rdkit.Contrib.SA_Score.sascorer",
        "numpy", "json", "math", "re"
    }

    # 금지된 키워드 (보안)
    FORBIDDEN_KEYWORDS = [
        "os.system", "subprocess", "eval(", "exec(",
        "open(", "__import__", "importlib",
        "socket", "urllib", "requests", "http",
        "shutil", "pathlib", "glob", "os.path",
        "os.remove", "os.rmdir", "os.unlink",
        "compile(", "getattr(", "setattr(",
    ]

    DEFAULT_TIMEOUT = 30  # seconds
    DOCKER_IMAGE = "adc-sandbox:latest"
    DOCKER_MEMORY_LIMIT = "512m"

    def __init__(self, mode: ExecutionMode = None):
        """
        초기화

        mode가 None이면 환경변수 SANDBOX_MODE 또는 기본값 사용
        """
        if mode:
            self.mode = mode
        else:
            env_mode = os.getenv("SANDBOX_MODE", "subprocess").lower()
            self.mode = ExecutionMode.DOCKER if env_mode == "docker" else ExecutionMode.SUBPROCESS

        logger.info(f"[sandbox] Initialized with mode: {self.mode.value}")

    def validate_code(self, code: str) -> tuple[bool, Optional[str]]:
        """
        코드 보안 검증

        Returns:
            (is_valid, error_message)
        """
        # 금지된 키워드 체크
        for keyword in self.FORBIDDEN_KEYWORDS:
            if keyword in code:
                return False, f"Forbidden keyword detected: {keyword}"

        # 위험한 import 체크
        import re
        import_pattern = r'(?:from\s+(\S+)\s+import|import\s+(\S+))'
        imports = re.findall(import_pattern, code)

        for imp in imports:
            module = imp[0] or imp[1]
            # 기본 모듈이나 허용된 모듈인지 확인
            root_module = module.split('.')[0]
            if root_module not in {'json', 'math', 're', 'rdkit', 'numpy'}:
                return False, f"Unauthorized import: {module}"

        return True, None

    async def execute(
        self,
        code: str,
        timeout: int = None,
        session_id: str = None
    ) -> ExecutionResult:
        """
        코드 실행

        Args:
            code: 실행할 Python 코드
            timeout: 타임아웃 (초)
            session_id: 세션 ID (로깅용)

        Returns:
            ExecutionResult
        """
        timeout = timeout or self.DEFAULT_TIMEOUT

        # 코드 검증
        is_valid, error = self.validate_code(code)
        if not is_valid:
            logger.warning(f"[sandbox] Code validation failed: {error}")
            return ExecutionResult(
                success=False,
                exit_code=-1,
                stdout="",
                stderr=f"Security validation failed: {error}",
                execution_time_ms=0,
                mode=self.mode.value
            )

        # 모드에 따라 실행
        if self.mode == ExecutionMode.DOCKER:
            return await self._execute_docker(code, timeout, session_id)
        else:
            return await self._execute_subprocess(code, timeout, session_id)

    @staticmethod
    def _get_preexec_fn():
        """리소스 제한 (Linux only) — Fat Container 보안 강화"""
        if platform.system() != "Linux":
            return None

        def _set_limits():
            import resource
            # 메모리 제한: 512MB
            mem_limit = 512 * 1024 * 1024
            resource.setrlimit(resource.RLIMIT_AS, (mem_limit, mem_limit))
            # CPU 시간 제한: 60초
            resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
            # 프로세스 생성 제한: 자식 프로세스 금지
            resource.setrlimit(resource.RLIMIT_NPROC, (0, 0))

        return _set_limits

    async def _execute_subprocess(
        self,
        code: str,
        timeout: int,
        session_id: str = None
    ) -> ExecutionResult:
        """subprocess로 직접 실행 (Fat Container 방식)"""
        import time
        start_time = time.time()

        # 임시 파일에 코드 저장
        with tempfile.NamedTemporaryFile(
            mode='w',
            suffix='.py',
            delete=False,
            encoding='utf-8'
        ) as f:
            f.write(code)
            script_path = f.name

        try:
            # 환경변수 전달 (PYTHONPATH, LD_LIBRARY_PATH 등 보장)
            env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1"}

            # sys.executable 사용 — Cloud Run에서 올바른 인터프리터 보장
            result = subprocess.run(
                [sys.executable, script_path],
                capture_output=True,
                text=True,
                timeout=timeout,
                env=env,
                preexec_fn=self._get_preexec_fn(),
            )

            execution_time_ms = int((time.time() - start_time) * 1000)

            if result.returncode != 0:
                logger.warning(
                    f"[sandbox] Script failed (exit={result.returncode}): "
                    f"stderr={result.stderr[:500]}"
                )

            return ExecutionResult(
                success=result.returncode == 0,
                exit_code=result.returncode,
                stdout=result.stdout,
                stderr=result.stderr,
                execution_time_ms=execution_time_ms,
                mode="subprocess"
            )

        except subprocess.TimeoutExpired:
            return ExecutionResult(
                success=False,
                exit_code=-1,
                stdout="",
                stderr=f"Execution timeout: {timeout}s exceeded",
                execution_time_ms=timeout * 1000,
                mode="subprocess"
            )
        except Exception as e:
            logger.exception(f"[sandbox] Subprocess error: {e}")
            return ExecutionResult(
                success=False,
                exit_code=-1,
                stdout="",
                stderr=str(e),
                execution_time_ms=0,
                mode="subprocess"
            )
        finally:
            # 임시 파일 삭제
            if os.path.exists(script_path):
                os.remove(script_path)

    async def _execute_docker(
        self,
        code: str,
        timeout: int,
        session_id: str = None
    ) -> ExecutionResult:
        """Docker 컨테이너에서 격리 실행 (프로덕션용)"""
        import time
        start_time = time.time()

        # 임시 디렉토리 생성
        with tempfile.TemporaryDirectory() as tmpdir:
            script_path = os.path.join(tmpdir, "script.py")

            with open(script_path, 'w', encoding='utf-8') as f:
                f.write(code)

            try:
                # Docker 실행 명령
                docker_cmd = [
                    "docker", "run",
                    "--rm",                              # 컨테이너 자동 삭제
                    "--network", "none",                 # 네트워크 차단
                    "--memory", self.DOCKER_MEMORY_LIMIT,# 메모리 제한
                    "--cpus", "1",                       # CPU 제한
                    "--read-only",                       # 읽기 전용 파일시스템
                    "--security-opt", "no-new-privileges",
                    "-v", f"{tmpdir}:/app:ro",           # 스크립트 마운트 (읽기전용)
                    "-w", "/app",
                    self.DOCKER_IMAGE,
                    "python", "/app/script.py"
                ]

                result = subprocess.run(
                    docker_cmd,
                    capture_output=True,
                    text=True,
                    timeout=timeout + 5  # Docker 오버헤드 고려
                )

                execution_time_ms = int((time.time() - start_time) * 1000)

                return ExecutionResult(
                    success=result.returncode == 0,
                    exit_code=result.returncode,
                    stdout=result.stdout,
                    stderr=result.stderr,
                    execution_time_ms=execution_time_ms,
                    mode="docker"
                )

            except subprocess.TimeoutExpired:
                # 타임아웃 시 컨테이너 강제 종료 시도
                subprocess.run(
                    ["docker", "kill", f"sandbox-{session_id}"],
                    capture_output=True
                )
                return ExecutionResult(
                    success=False,
                    exit_code=-1,
                    stdout="",
                    stderr=f"Docker execution timeout: {timeout}s exceeded",
                    execution_time_ms=timeout * 1000,
                    mode="docker"
                )
            except FileNotFoundError:
                # Docker가 설치되지 않은 경우 subprocess로 폴백
                logger.warning("[sandbox] Docker not found, falling back to subprocess")
                return await self._execute_subprocess(code, timeout, session_id)
            except Exception as e:
                logger.exception(f"[sandbox] Docker error: {e}")
                return ExecutionResult(
                    success=False,
                    exit_code=-1,
                    stdout="",
                    stderr=str(e),
                    execution_time_ms=0,
                    mode="docker"
                )

    async def check_docker_available(self) -> bool:
        """Docker 사용 가능 여부 확인"""
        try:
            result = subprocess.run(
                ["docker", "version"],
                capture_output=True,
                timeout=5
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False

    async def check_image_exists(self) -> bool:
        """Sandbox 이미지 존재 여부 확인"""
        try:
            result = subprocess.run(
                ["docker", "images", "-q", self.DOCKER_IMAGE],
                capture_output=True,
                text=True,
                timeout=5
            )
            return bool(result.stdout.strip())
        except Exception:
            return False


    async def verify_rdkit_available(self) -> tuple[bool, str]:
        """
        RDKit 설치 상태 검증 (startup health check용)

        Returns:
            (is_available, message)
        """
        test_code = """
import json
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    mol = Chem.MolFromSmiles('CCO')
    assert mol is not None, "MolFromSmiles returned None"
    mw = Descriptors.MolWt(mol)
    result = {"ok": True, "rdkit_version": Chem.rdBase.rdkitVersion, "test_mw": round(mw, 2)}
except Exception as e:
    result = {"ok": False, "error": str(e)}
print(json.dumps(result))
"""
        exec_result = await self.execute(test_code, timeout=15)
        if exec_result.success:
            try:
                data = json.loads(exec_result.stdout)
                if data.get("ok"):
                    return True, f"RDKit {data.get('rdkit_version', '?')} — test MW(CCO)={data.get('test_mw')}"
                else:
                    return False, f"RDKit test failed: {data.get('error')}"
            except json.JSONDecodeError:
                return False, f"Unexpected output: {exec_result.stdout[:200]}"
        else:
            return False, f"RDKit verification script failed: {exec_result.stderr[:300]}"


# 싱글톤 인스턴스
_sandbox_executor: Optional[SandboxExecutor] = None


def get_sandbox_executor() -> SandboxExecutor:
    """Sandbox Executor 싱글톤 반환"""
    global _sandbox_executor
    if _sandbox_executor is None:
        _sandbox_executor = SandboxExecutor()
    return _sandbox_executor
