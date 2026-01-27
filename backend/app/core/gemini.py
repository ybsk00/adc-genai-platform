"""
Gemini Client Module
LangChain + Google Gemini 연동을 위한 클라이언트 모듈
"""
from langchain_google_genai import ChatGoogleGenerativeAI
from app.core.config import settings


def get_gemini_model(
    temperature: float = 0,
    model_name: str | None = None
) -> ChatGoogleGenerativeAI:
    """
    Gemini 2.5 Flash 모델 인스턴스 반환
    
    Args:
        temperature: 생성 온도 (0~1)
        model_name: 모델 이름 (기본값: settings.GEMINI_MODEL_ID)
    
    Returns:
        ChatGoogleGenerativeAI: Gemini 채팅 모델
    """
    return ChatGoogleGenerativeAI(
        model=model_name or settings.GEMINI_MODEL_ID,
        temperature=temperature,
        google_api_key=settings.GOOGLE_API_KEY,
        convert_system_message_to_human=True  # Gemini는 system role 미지원시 변환
    )


# Singleton instance for common use
# gemini_flash = get_gemini_model(temperature=0)  # Removed to prevent startup crash if API key is missing

