from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from app.api import auth, jobs, admin, library, payment, scheduler, knowledge_base
from app.core.config import settings

app = FastAPI(
    title="ADC-GenAI Platform API",
    description="AI-driven ADC analysis platform for bio-researchers",
    version="1.0.0"
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


@app.get("/health")
async def health_check():
    """헬스 체크 엔드포인트"""
    return {"status": "healthy", "version": "1.0.0"}
