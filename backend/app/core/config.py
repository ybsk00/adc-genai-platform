from pydantic_settings import BaseSettings
from typing import List


class Settings(BaseSettings):
    """애플리케이션 설정"""
    
    # App
    APP_NAME: str = "ADC-GenAI Platform"
    DEBUG: bool = False
    
    # CORS
    CORS_ORIGINS: List[str] = [
        "http://localhost:5173", 
        "http://localhost:3000",
        "https://astraforge-adc-platform.web.app",
        "https://astraforge-adc-platform.firebaseapp.com"
    ]

    
    # Supabase
    SUPABASE_URL: str = ""
    SUPABASE_SERVICE_KEY: str = ""
    
    # OpenAI
    OPENAI_API_KEY: str = ""

    # [Model Strategy] 가성비 전략
    # 리포트 작성, 심층 구조 분석용 (Brain)
    SMART_LLM: str = "gpt-4o"
    
    # 단순 검색, 데이터 정리, 분류용 (Worker) - 10배 저렴
    FAST_LLM: str = "gpt-4o-mini"

    
    # Perplexity
    PERPLEXITY_API_KEY: str = ""
    
    # BioNeMo
    BIONEMO_API_KEY: str = ""
    
    # Tavily
    TAVILY_API_KEY: str = ""
    
    # Google (Gemini)
    GOOGLE_API_KEY: str = ""
    GEMINI_MODEL_ID: str = "gemini-2.0-flash"

    # NCBI (PubMed)
    NCBI_EMAIL: str = "your-email@example.com"
    NCBI_API_KEY: str = ""
    NCBI_TOOL: str = "adc-platform"
    
    # JWT
    JWT_SECRET_KEY: str = "your-secret-key-change-in-production"
    JWT_ALGORITHM: str = "HS256"
    ACCESS_TOKEN_EXPIRE_MINUTES: int = 30
    
    # Proxy
    PROXY_ENABLED: bool = False
    PROXY_HOST: str = ""
    PROXY_PORTS: str = ""
    PROXY_USERNAME: str = ""
    PROXY_PASSWORD: str = ""
    PROXY_LOCATION: str = "us"
    PROXY_ROTATING: bool = True
    
    class Config:
        env_file = ".env"
        case_sensitive = True
        extra = "ignore"


settings = Settings()
