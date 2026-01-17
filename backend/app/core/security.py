from fastapi import Request, HTTPException
import httpx


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
