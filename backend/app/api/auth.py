from fastapi import APIRouter

router = APIRouter()


@router.post("/lead-magnet")
async def collect_lead_email(email: str):
    """랜딩페이지 이메일 수집 (Golden Set PDF 다운로드 링크 발송)"""
    # TODO: CRM DB에 이메일 저장 및 자동 이메일 발송
    return {"status": "sent", "msg": "Check inbox for Golden Set download link"}


@router.get("/profile")
async def get_user_profile():
    """사용자 정보 및 크레딧 조회"""
    # TODO: JWT에서 사용자 ID 추출 후 DB 조회
    return {
        "id": "u_1",
        "email": "user@example.com",
        "credits": 50,
        "plan": "free"
    }
