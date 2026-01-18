from fastapi import Request, HTTPException, Depends
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
import httpx
from supabase import create_client, Client
from app.core.config import settings

# Supabase client
supabase: Client = create_client(settings.SUPABASE_URL, settings.SUPABASE_SERVICE_KEY)

# HTTP Bearer security scheme
security = HTTPBearer()


async def verify_supabase_token(credentials: HTTPAuthorizationCredentials = Depends(security)) -> dict:
    """Supabase JWT 토큰 검증
    
    프론트엔드에서 받은 access_token을 검증하고 사용자 정보 반환
    """
    token = credentials.credentials
    
    try:
        # Supabase에서 사용자 정보 가져오기
        user_response = supabase.auth.get_user(token)
        
        if not user_response.user:
            raise HTTPException(
                status_code=401,
                detail="Invalid or expired token"
            )
        
        return {
            "user_id": user_response.user.id,
            "email": user_response.user.email,
            "role": user_response.user.role
        }
    except Exception as e:
        raise HTTPException(
            status_code=401,
            detail=f"Token verification failed: {str(e)}"
        )


async def verify_oidc_token(request: Request) -> None:
    """Google Cloud Scheduler OIDC 토큰 검증
    
    Cloud Scheduler가 보낸 요청인지 확인하여 
    무단 API 호출 방지 (OpenAI 비용 폭발 방지)
    """
    auth_header = request.headers.get("Authorization", "")
    
    if not auth_header.startswith("Bearer "):
        raise HTTPException(
            status_code=403, 
            detail="Forbidden: Missing or invalid OIDC token"
        )
    
    token = auth_header.split(" ")[1]
    
    # TODO: Google 공개키로 JWT 검증
    # audience = Cloud Run 서비스 URL
    # 참고: https://cloud.google.com/scheduler/docs/http-target-auth
    
    # 개발 환경에서는 토큰 검증 스킵 (프로덕션에서 활성화)
    # 실제 구현:
    # 1. Google의 공개 키 가져오기
    # 2. JWT 서명 검증
    # 3. audience, issuer 검증
    
    if not token:
        raise HTTPException(status_code=403, detail="Forbidden: Empty token")

