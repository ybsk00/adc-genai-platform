import os
# Fix for gRPC fork support (must be set before any grpc import)
os.environ["GRPC_ENABLE_FORK_SUPPORT"] = "false"
os.environ["GRPC_POLL_STRATEGY"] = "epoll1"
os.environ["GRPC_TYPE_CHECK"] = "0"

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager
import logging

# 로깅 설정 (INFO 레벨)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

from app.api import auth, jobs, admin, library, payment, scheduler, knowledge_base, design_runs, assay_results, design, uniprot, system, report
from app.core.config import settings
from app.services.scheduler_engine import scheduler_engine
from app.services.sandbox_executor import get_sandbox_executor

@asynccontextmanager
async def lifespan(app: FastAPI):
    # 앱 시작 시 스케줄러 가동
    scheduler_engine.start()

    # Sandbox + RDKit 상태 검증 (Fat Container Health Check)
    try:
        sandbox = get_sandbox_executor()
        ok, msg = await sandbox.verify_rdkit_available()
        if ok:
            logger.info(f"[startup] Sandbox RDKit verified: {msg}")
        else:
            logger.critical(f"[startup] Sandbox RDKit UNAVAILABLE: {msg}")
    except Exception as e:
        logger.critical(f"[startup] Sandbox verification failed: {e}")

    yield
    # 앱 종료 시 스케줄러 중지
    scheduler_engine.shutdown()

app = FastAPI(
    title="ADC-GenAI Platform API",
    description="AI-driven ADC analysis platform for bio-researchers",
    version="1.0.0",
    lifespan=lifespan
)

# CORS 설정
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# 라우터 등록
app.include_router(auth.router, prefix="/api/auth", tags=["auth"])
app.include_router(jobs.router, prefix="/api/jobs", tags=["jobs"])
app.include_router(admin.router, prefix="/api/admin", tags=["admin"])
app.include_router(library.router, prefix="/api/library", tags=["library"])
app.include_router(payment.router, prefix="/api/payment", tags=["payment"])
app.include_router(payment.router, prefix="/api/payment", tags=["payment"])
app.include_router(scheduler.router, prefix="/api/scheduler", tags=["scheduler"])
app.include_router(knowledge_base.router, prefix="/api/knowledge-base", tags=["knowledge-base"])
app.include_router(design_runs.router, prefix="/api/design-runs", tags=["design-runs"])
app.include_router(assay_results.router, prefix="/api/assay-results", tags=["assay-results"])
app.include_router(design.router, tags=["design-engine"])  # prefix already in router
app.include_router(uniprot.router, prefix="/api/uniprot", tags=["uniprot"])
app.include_router(system.router, prefix="/api/system", tags=["system"])
app.include_router(report.router, prefix="/api/report", tags=["report"])


@app.get("/health")
async def health_check():
    """헬스 체크 엔드포인트"""
    return {"status": "healthy", "version": "1.0.0"}
