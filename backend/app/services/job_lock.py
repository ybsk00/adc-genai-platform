"""
Job Lock Service
중복 작업 실행 방지를 위한 DB 기반 잠금 메커니즘
"""
import logging
from datetime import datetime, timedelta
from typing import Optional

from app.core.supabase import supabase

logger = logging.getLogger(__name__)


class JobLock:
    """DB 기반 작업 잠금 메커니즘"""
    
    # 잠금 만료 시간 (시간 단위) - 비정상 종료 시 자동 해제
    LOCK_EXPIRY_HOURS = 6
    
    async def acquire(self, job_type: str) -> bool:
        """
        작업 잠금 획득 시도
        - 실행 중인 같은 타입의 작업이 없으면 True 반환
        - 이미 실행 중이면 False 반환
        """
        try:
            # 1. 오래된 잠금 정리 (비정상 종료 대비)
            expiry_time = datetime.utcnow() - timedelta(hours=self.LOCK_EXPIRY_HOURS)
            supabase.table("job_locks").delete().lt("acquired_at", expiry_time.isoformat()).execute()
            
            # 2. 현재 잠금 확인
            existing = supabase.table("job_locks")\
                .select("id, acquired_at")\
                .eq("job_type", job_type)\
                .execute()
            
            if existing.data:
                lock = existing.data[0]
                logger.warning(f"⚠️ Job '{job_type}' is already running since {lock['acquired_at']}")
                return False
            
            # 3. 새 잠금 생성
            supabase.table("job_locks").insert({
                "job_type": job_type,
                "acquired_at": datetime.utcnow().isoformat()
            }).execute()
            
            logger.info(f"🔒 Lock acquired for job: {job_type}")
            return True
            
        except Exception as e:
            logger.error(f"Lock acquire error: {e}")
            # 잠금 실패 시에도 작업 실행 허용 (fail-open)
            return True
    
    async def release(self, job_type: str):
        """작업 잠금 해제"""
        try:
            supabase.table("job_locks").delete().eq("job_type", job_type).execute()
            logger.info(f"🔓 Lock released for job: {job_type}")
        except Exception as e:
            logger.error(f"Lock release error: {e}")
    
    async def is_locked(self, job_type: str) -> bool:
        """잠금 상태 확인"""
        try:
            existing = supabase.table("job_locks")\
                .select("id")\
                .eq("job_type", job_type)\
                .execute()
            return bool(existing.data)
        except Exception:
            return False

    # --- Component Locking (Async-Precompute) ---

    async def acquire_component_lock(self, component_id: str, worker_id: str, estimated_wait: int = 60) -> bool:
        """
        컴포넌트(분자) 단위 잠금 획득
        - component_catalog 테이블의 lock_status 업데이트
        """
        try:
            # 1. 현재 상태 확인
            res = supabase.table("component_catalog")\
                .select("lock_status, lock_holder")\
                .eq("id", component_id)\
                .execute()
            
            if not res.data:
                return False # 컴포넌트가 없음
            
            item = res.data[0]
            if item.get("lock_status") == "computing":
                # 이미 계산 중
                logger.info(f"🔒 Component {component_id} is already being computed by {item.get('lock_holder')}")
                return False
            
            # 2. 잠금 시도 (Update)
            update_data = {
                "lock_status": "computing",
                "lock_holder": worker_id,
                "estimated_wait_time": estimated_wait,
                "updated_at": datetime.utcnow().isoformat()
            }
            
            supabase.table("component_catalog")\
                .update(update_data)\
                .eq("id", component_id)\
                .execute()
                
            return True
            
        except Exception as e:
            logger.error(f"Component lock error: {e}")
            return False

    async def release_component_lock(self, component_id: str):
        """컴포넌트 잠금 해제"""
        try:
            update_data = {
                "lock_status": "available",
                "lock_holder": None,
                "estimated_wait_time": 0,
                "updated_at": datetime.utcnow().isoformat()
            }
            
            supabase.table("component_catalog")\
                .update(update_data)\
                .eq("id", component_id)\
                .execute()
                
        except Exception as e:
            logger.error(f"Component unlock error: {e}")


# 싱글톤 인스턴스
job_lock = JobLock()
