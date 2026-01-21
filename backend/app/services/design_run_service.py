"""
Design Run Service
Manages design runs with frozen_params for reproducibility (Eng-Fit v0.2)
"""
import logging
from datetime import datetime
from typing import Optional, Dict, Any
from uuid import uuid4

from app.core.supabase import supabase

logger = logging.getLogger(__name__)


# Default scoring parameters (Eng-Fit v0.2)
# These can be fetched from system_config or a dedicated scoring_params table
DEFAULT_SCORING_PARAMS = {
    "version": "scoring_v0.2",
    "weights": {
        "efficacy": 0.35,
        "safety": 0.25,
        "manufacturability": 0.20,
        "ip_freedom": 0.10,
        "cost": 0.10
    },
    "thresholds": {
        "agg_threshold": 5.0,  # % aggregation
        "binding_affinity_min": 1.0,  # nM
        "dar_optimal_range": [3.5, 4.5]
    }
}


class DesignRunService:
    """
    Design Run 관리 서비스
    - frozen_params 스냅샷 저장
    - 재현성 보장 (Audit 가능)
    """

    async def get_current_scoring_params(self) -> Dict[str, Any]:
        """
        현재 활성 스코링 파라미터 조회
        system_config 테이블 또는 기본값 반환
        """
        try:
            res = supabase.table("system_config")\
                .select("value")\
                .eq("key", "SCORING_PARAMS")\
                .execute()
            
            if res.data:
                return res.data[0]["value"]
            return DEFAULT_SCORING_PARAMS
        except Exception as e:
            logger.warning(f"Failed to fetch scoring params, using default: {e}")
            return DEFAULT_SCORING_PARAMS

    async def create_design_run(
        self, 
        project_id: str, 
        user_id: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        새 Design Run 생성 (frozen_params 스냅샷 포함)
        """
        try:
            # 1. 현재 스코링 파라미터 스냅샷
            frozen_params = await self.get_current_scoring_params()
            frozen_params["frozen_at"] = datetime.utcnow().isoformat()
            
            # 2. Design Run 레코드 생성
            run_id = str(uuid4())
            new_run = {
                "id": run_id,
                "project_id": project_id,
                "status": "draft",
                "frozen_params": frozen_params,
                "created_at": datetime.utcnow().isoformat(),
                "updated_at": datetime.utcnow().isoformat()
            }
            
            res = supabase.table("design_runs").insert(new_run).execute()
            
            if res.data:
                logger.info(f"✅ Created design run {run_id} with frozen_params v{frozen_params.get('version')}")
                return res.data[0]
            else:
                raise Exception("Failed to insert design_run")
                
        except Exception as e:
            logger.error(f"Create design run error: {e}")
            raise

    async def update_run_status(
        self, 
        run_id: str, 
        status: str
    ):
        """Design Run 상태 업데이트"""
        try:
            update_data = {
                "status": status,
                "updated_at": datetime.utcnow().isoformat()
            }
            supabase.table("design_runs").update(update_data).eq("id", run_id).execute()
            logger.info(f"🔄 Design run {run_id} status -> {status}")
        except Exception as e:
            logger.error(f"Update run status error: {e}")

    async def get_run(self, run_id: str) -> Optional[Dict[str, Any]]:
        """Design Run 조회"""
        try:
            res = supabase.table("design_runs").select("*").eq("id", run_id).execute()
            return res.data[0] if res.data else None
        except Exception as e:
            logger.error(f"Get run error: {e}")
            return None


# Singleton
design_run_service = DesignRunService()
