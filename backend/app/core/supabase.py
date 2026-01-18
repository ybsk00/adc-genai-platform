"""
Supabase Client
"""
from supabase import create_client, Client
from app.core.config import settings

def get_supabase_client() -> Client:
    """Supabase 클라이언트 반환"""
    url = settings.SUPABASE_URL
    key = settings.SUPABASE_SERVICE_KEY
    
    if not url or not key:
        print("Warning: SUPABASE_URL or SUPABASE_SERVICE_KEY is missing.")
        # 개발 환경에서 에러 방지를 위해 더미 클라이언트나 에러 처리 필요
        # 여기서는 일단 진행
    
    return create_client(url, key)

supabase = get_supabase_client()
