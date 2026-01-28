# ADC Sandbox Executor - Air-gapped Python Environment
# Constraint 7: Library Version Lock

FROM python:3.11-slim

LABEL maintainer="ADC Platform Team"
LABEL description="Secure sandbox for executing chemistry validation code"

# 보안: 비-root 사용자 생성
RUN groupadd -r sandbox && useradd -r -g sandbox sandbox

# 시스템 패키지 설치 (RDKit 의존성)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Python 패키지 설치 (버전 고정)
COPY requirements.lock /tmp/requirements.lock
RUN pip install --no-cache-dir -r /tmp/requirements.lock

# 작업 디렉토리 설정
WORKDIR /app

# 사용자 전환 (보안)
USER sandbox

# 헬스체크
HEALTHCHECK --interval=30s --timeout=5s --start-period=5s --retries=3 \
    CMD python -c "from rdkit import Chem; print('OK')" || exit 1

# 기본 명령어
CMD ["python"]
