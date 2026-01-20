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
    && rm -rf /var/lib/apt/lists/*

# 작업 디렉토리
WORKDIR /app

# 의존성 설치
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 앱 코드 복사
COPY app/ ./app/

# 환경 변수
ENV PORT=8080
ENV PYTHONUNBUFFERED=1

# 실행
# 시작 스크립트 복사 및 실행 권한 부여
COPY start.sh .
RUN chmod +x start.sh

# 실행
CMD ["./start.sh"]
