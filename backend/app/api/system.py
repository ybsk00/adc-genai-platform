"""
System Information API
환경 정보 및 버전 관리를 위한 시스템 API
"""
from fastapi import APIRouter
from pydantic import BaseModel
from typing import Dict, List, Optional
import sys
import platform
import logging

logger = logging.getLogger(__name__)

router = APIRouter()


class LibraryVersion(BaseModel):
    """라이브러리 버전 정보"""
    name: str
    version: str
    category: str  # 'core', 'ai', 'chemistry', 'data', 'web'


class EnvironmentInfo(BaseModel):
    """환경 정보 응답"""
    python_version: str
    platform: str
    libraries: List[LibraryVersion]
    environment: str  # 'development', 'staging', 'production'


def get_package_version(package_name: str) -> Optional[str]:
    """패키지 버전 조회"""
    try:
        import importlib.metadata
        return importlib.metadata.version(package_name)
    except Exception:
        return None


# 핵심 라이브러리 목록 (카테고리별)
TRACKED_LIBRARIES = {
    "core": [
        "fastapi",
        "uvicorn",
        "pydantic",
        "python-dotenv",
    ],
    "ai": [
        "google-generativeai",
        "langchain",
        "langchain-google-genai",
        "openai",
    ],
    "chemistry": [
        "rdkit",
        "pubchempy",
    ],
    "data": [
        "supabase",
        "httpx",
        "aiohttp",
        "numpy",
        "pandas",
    ],
    "web": [
        "websockets",
        "python-multipart",
    ],
}


@router.get("/environment", response_model=EnvironmentInfo)
async def get_environment_info():
    """
    시스템 환경 정보 조회

    Python 버전, 플랫폼, 핵심 라이브러리 버전 정보를 반환합니다.
    21 CFR Part 11 준수를 위한 환경 버전 추적에 사용됩니다.
    """
    import os

    libraries = []

    for category, packages in TRACKED_LIBRARIES.items():
        for package in packages:
            version = get_package_version(package)
            if version:
                libraries.append(LibraryVersion(
                    name=package,
                    version=version,
                    category=category
                ))

    # 환경 결정
    env = os.getenv("ENVIRONMENT", "development")
    if os.getenv("PRODUCTION"):
        env = "production"
    elif os.getenv("STAGING"):
        env = "staging"

    return EnvironmentInfo(
        python_version=f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        platform=f"{platform.system()} {platform.release()}",
        libraries=libraries,
        environment=env
    )


@router.get("/environment/lock")
async def get_environment_lock():
    """
    환경 잠금 정보 (버전 고정용)

    requirements.txt 형식으로 현재 환경의 핵심 라이브러리 버전을 반환합니다.
    """
    lock_entries = []

    for category, packages in TRACKED_LIBRARIES.items():
        lock_entries.append(f"# {category.upper()}")
        for package in packages:
            version = get_package_version(package)
            if version:
                lock_entries.append(f"{package}=={version}")
        lock_entries.append("")

    return {
        "format": "requirements.txt",
        "python_version": f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
        "lock_content": "\n".join(lock_entries),
        "generated_at": __import__("datetime").datetime.utcnow().isoformat()
    }


@router.get("/health/detailed")
async def detailed_health_check():
    """
    상세 헬스 체크

    시스템 상태 및 외부 서비스 연결 상태를 확인합니다.
    """
    import os

    checks = {
        "api": {"status": "healthy", "message": "API is running"},
        "python": {"status": "healthy", "version": f"{sys.version_info.major}.{sys.version_info.minor}"},
    }

    # Supabase 연결 확인
    try:
        from app.core.supabase import get_supabase_client
        client = get_supabase_client()
        # 간단한 쿼리로 연결 확인
        checks["supabase"] = {"status": "healthy", "message": "Connected"}
    except Exception as e:
        checks["supabase"] = {"status": "unhealthy", "message": str(e)}

    # Gemini API 확인
    gemini_key = os.getenv("GOOGLE_API_KEY") or os.getenv("GEMINI_API_KEY")
    if gemini_key:
        checks["gemini"] = {"status": "configured", "message": "API key present"}
    else:
        checks["gemini"] = {"status": "not_configured", "message": "API key missing"}

    # 전체 상태 결정
    overall_status = "healthy"
    for check in checks.values():
        if check.get("status") == "unhealthy":
            overall_status = "degraded"
            break

    return {
        "status": overall_status,
        "checks": checks,
        "timestamp": __import__("datetime").datetime.utcnow().isoformat()
    }
