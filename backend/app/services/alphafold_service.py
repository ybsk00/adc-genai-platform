"""
AlphaFold 3 Integration Service
Premium 사용자 전용 3D 구조 예측 서비스

Constraint: GPU Quota Management
- Monthly quota: 50 predictions per premium user
- Daily quota: 5 predictions
"""
import asyncio
import logging
from typing import Optional, Dict, Any
from datetime import datetime, timedelta
from enum import Enum

from app.core.supabase import get_supabase_client

logger = logging.getLogger(__name__)


class AlphaFoldJobStatus(Enum):
    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"


class QuotaExceededError(Exception):
    """사용자 할당량 초과 예외"""
    pass


class AlphaFoldService:
    """
    AlphaFold 3 서비스

    Phase 3.5: Premium 사용자 전용 3D 구조 예측
    - Google Cloud AlphaFold API 연동
    - 사용자별 GPU 쿼터 관리
    - 비동기 작업 처리
    """

    def __init__(self):
        self.supabase = get_supabase_client()
        # AlphaFold API 설정 (향후 구현)
        self.api_endpoint = None
        self.api_key = None

    async def check_quota(self, user_id: str) -> Dict[str, Any]:
        """
        사용자 쿼터 확인

        Args:
            user_id: 사용자 ID

        Returns:
            {
                "can_submit": bool,
                "monthly_used": int,
                "monthly_limit": int,
                "daily_used": int,
                "daily_limit": int,
                "reset_at": datetime
            }
        """
        try:
            # 쿼터 테이블 조회
            result = await self.supabase.table("alphafold_user_quota").select("*").eq(
                "user_id", user_id
            ).single().execute()

            if not result.data:
                # 쿼터 레코드가 없으면 생성
                await self._create_quota_record(user_id)
                return {
                    "can_submit": True,
                    "monthly_used": 0,
                    "monthly_limit": 50,
                    "daily_used": 0,
                    "daily_limit": 5,
                    "reset_at": self._get_next_month_reset()
                }

            quota = result.data

            # 월 리셋 확인
            if self._should_reset_monthly(quota.get("last_reset_at")):
                await self._reset_monthly_quota(user_id)
                quota["monthly_used"] = 0

            # 일일 리셋 확인
            if self._should_reset_daily(quota.get("last_daily_reset")):
                await self._reset_daily_quota(user_id)
                quota["daily_used"] = 0

            monthly_used = quota.get("monthly_used", 0)
            monthly_limit = quota.get("monthly_limit", 50)
            daily_used = quota.get("daily_used", 0)
            daily_limit = quota.get("daily_limit", 5)

            can_submit = monthly_used < monthly_limit and daily_used < daily_limit

            return {
                "can_submit": can_submit,
                "monthly_used": monthly_used,
                "monthly_limit": monthly_limit,
                "daily_used": daily_used,
                "daily_limit": daily_limit,
                "reset_at": self._get_next_month_reset()
            }

        except Exception as e:
            logger.error(f"Failed to check quota: {e}")
            raise

    async def submit_prediction(
        self,
        user_id: str,
        session_id: str,
        sequence: str,
        structure_type: str = "antibody"
    ) -> str:
        """
        AlphaFold 예측 작업 제출

        Args:
            user_id: 사용자 ID
            session_id: 설계 세션 ID
            sequence: 아미노산 서열
            structure_type: 구조 타입 (antibody, payload, linker)

        Returns:
            job_id: 생성된 작업 ID

        Raises:
            QuotaExceededError: 쿼터 초과 시
        """
        # 쿼터 확인
        quota = await self.check_quota(user_id)
        if not quota["can_submit"]:
            raise QuotaExceededError(
                f"Quota exceeded. Monthly: {quota['monthly_used']}/{quota['monthly_limit']}, "
                f"Daily: {quota['daily_used']}/{quota['daily_limit']}"
            )

        try:
            # 작업 생성
            result = await self.supabase.table("alphafold_jobs").insert({
                "session_id": session_id,
                "user_id": user_id,
                "input_sequence": sequence,
                "structure_type": structure_type,
                "status": AlphaFoldJobStatus.PENDING.value
            }).execute()

            if not result.data:
                raise ValueError("Failed to create AlphaFold job")

            job_id = result.data[0]["id"]

            # 쿼터 증가
            await self._increment_quota(user_id)

            # 비동기 작업 큐에 추가 (실제 구현 시 Celery 등 사용)
            asyncio.create_task(self._process_prediction(job_id, sequence, structure_type))

            logger.info(f"AlphaFold job submitted: {job_id}")
            return job_id

        except QuotaExceededError:
            raise
        except Exception as e:
            logger.error(f"Failed to submit prediction: {e}")
            raise

    async def get_job_status(self, job_id: str, user_id: str) -> Dict[str, Any]:
        """
        작업 상태 조회

        Args:
            job_id: 작업 ID
            user_id: 사용자 ID (권한 확인용)

        Returns:
            작업 상태 정보
        """
        result = await self.supabase.table("alphafold_jobs").select("*").eq(
            "id", job_id
        ).eq(
            "user_id", user_id
        ).single().execute()

        if not result.data:
            raise ValueError("Job not found or access denied")

        return result.data

    async def get_prediction_result(self, job_id: str, user_id: str) -> Optional[Dict]:
        """
        예측 결과 조회

        Args:
            job_id: 작업 ID
            user_id: 사용자 ID

        Returns:
            PDB 구조 및 신뢰도 점수
        """
        job = await self.get_job_status(job_id, user_id)

        if job["status"] != AlphaFoldJobStatus.COMPLETED.value:
            return None

        return {
            "pdb_content": job.get("output_pdb"),
            "confidence_score": job.get("plddt_score"),
            "predicted_tm_score": job.get("ptm_score"),
            "model_version": job.get("model_version")
        }

    async def _process_prediction(
        self,
        job_id: str,
        sequence: str,
        structure_type: str
    ):
        """
        실제 예측 처리 (비동기)

        TODO: Google Cloud AlphaFold API 연동
        현재는 플레이스홀더
        """
        try:
            # 상태 업데이트: RUNNING
            await self.supabase.table("alphafold_jobs").update({
                "status": AlphaFoldJobStatus.RUNNING.value,
                "started_at": datetime.utcnow().isoformat()
            }).eq("id", job_id).execute()

            # TODO: 실제 AlphaFold API 호출
            # 현재는 시뮬레이션
            await asyncio.sleep(5)  # 시뮬레이션 대기

            # 상태 업데이트: COMPLETED (플레이스홀더)
            await self.supabase.table("alphafold_jobs").update({
                "status": AlphaFoldJobStatus.COMPLETED.value,
                "completed_at": datetime.utcnow().isoformat(),
                "output_pdb": "PLACEHOLDER_PDB_CONTENT",
                "plddt_score": 85.0,  # 플레이스홀더
                "ptm_score": 0.9,
                "model_version": "alphafold3-preview"
            }).eq("id", job_id).execute()

            logger.info(f"AlphaFold job completed: {job_id}")

        except Exception as e:
            logger.error(f"AlphaFold prediction failed: {e}")

            await self.supabase.table("alphafold_jobs").update({
                "status": AlphaFoldJobStatus.FAILED.value,
                "error_message": str(e)
            }).eq("id", job_id).execute()

    async def _create_quota_record(self, user_id: str):
        """새 쿼터 레코드 생성"""
        await self.supabase.table("alphafold_user_quota").insert({
            "user_id": user_id,
            "monthly_limit": 50,
            "daily_limit": 5,
            "monthly_used": 0,
            "daily_used": 0,
            "last_reset_at": datetime.utcnow().isoformat(),
            "last_daily_reset": datetime.utcnow().isoformat()
        }).execute()

    async def _increment_quota(self, user_id: str):
        """쿼터 사용량 증가"""
        # 현재 값 조회
        result = await self.supabase.table("alphafold_user_quota").select(
            "monthly_used, daily_used"
        ).eq("user_id", user_id).single().execute()

        if result.data:
            await self.supabase.table("alphafold_user_quota").update({
                "monthly_used": result.data["monthly_used"] + 1,
                "daily_used": result.data["daily_used"] + 1
            }).eq("user_id", user_id).execute()

    async def _reset_monthly_quota(self, user_id: str):
        """월간 쿼터 리셋"""
        await self.supabase.table("alphafold_user_quota").update({
            "monthly_used": 0,
            "last_reset_at": datetime.utcnow().isoformat()
        }).eq("user_id", user_id).execute()

    async def _reset_daily_quota(self, user_id: str):
        """일일 쿼터 리셋"""
        await self.supabase.table("alphafold_user_quota").update({
            "daily_used": 0,
            "last_daily_reset": datetime.utcnow().isoformat()
        }).eq("user_id", user_id).execute()

    def _should_reset_monthly(self, last_reset: Optional[str]) -> bool:
        """월간 리셋 필요 여부"""
        if not last_reset:
            return True

        last_reset_dt = datetime.fromisoformat(last_reset.replace("Z", "+00:00"))
        now = datetime.utcnow()

        # 다른 월이면 리셋
        return last_reset_dt.month != now.month or last_reset_dt.year != now.year

    def _should_reset_daily(self, last_reset: Optional[str]) -> bool:
        """일일 리셋 필요 여부"""
        if not last_reset:
            return True

        last_reset_dt = datetime.fromisoformat(last_reset.replace("Z", "+00:00"))
        now = datetime.utcnow()

        # 다른 날이면 리셋
        return last_reset_dt.date() != now.date()

    def _get_next_month_reset(self) -> datetime:
        """다음 월 리셋 시간"""
        now = datetime.utcnow()
        if now.month == 12:
            return datetime(now.year + 1, 1, 1)
        return datetime(now.year, now.month + 1, 1)
