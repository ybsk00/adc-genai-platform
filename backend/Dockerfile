# ADC-GenAI Platform Backend
# Cloud Run용 Docker 이미지

FROM python:3.12-slim-bookworm

# WeasyPrint 의존성 설치
RUN apt-get update && apt-get install -y \
    build-essential \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libgdk-pixbuf2.0-0 \
    libffi-dev \
    shared-mime-info \
    fonts-noto-cjk \
    curl \
    libxrender1 \
    libxext6 \
    libgl1 \
    && rm -rf /var/lib/apt/lists/*

# 작업 디렉토리
WORKDIR /app

# 의존성 설치
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
RUN playwright install --with-deps chromium

# RDKit 설치 검증 (빌드 시점에 확인 — 실패 시 Docker build 중단)
RUN python -c "from rdkit import Chem; mol=Chem.MolFromSmiles('CCO'); assert mol is not None; print('RDKit OK:', Chem.rdBase.rdkitVersion)"

# 앱 코드 복사
COPY app/ ./app/
COPY run_crawler.py .

# 환경 변수
ENV PORT=8080
ENV PYTHONUNBUFFERED=1
ENV GRPC_ENABLE_FORK_SUPPORT=false

# 실행
# 시작 스크립트 복사 및 실행 권한 부여
COPY start.sh .
RUN chmod +x start.sh

# 실행
CMD ["./start.sh"]
