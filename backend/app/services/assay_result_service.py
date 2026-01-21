"""
Assay Result Service
Manages assay results with flexible success criteria (Closed-loop)
"""
import logging
from datetime import datetime
from typing import Optional, Dict, Any, List
from uuid import uuid4

from app.core.supabase import supabase

logger = logging.getLogger(__name__)


class AssayResultService:
    """
    Assay Result 관리 서비스
    - acceptance_criteria와 함께 결과 저장
    - Rule Confidence 계산을 위한 유연한 성공 판정
    """

    async def save_assay_result(
        self,
        run_id: str,
        molecule_id: str,
        assay_type: str,
        raw_data: Dict[str, Any],
        acceptance_criteria: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Assay 결과 저장 (acceptance_criteria 포함)
        
        Args:
            run_id: Design Run ID
            molecule_id: 분자 ID (golden_set_library 또는 component_catalog)
            assay_type: 분석 유형 (aggregation, binding, cytotoxicity 등)
            raw_data: 실험 데이터 (예: {"aggregation_percent": 3.2, "ph": 7.4})
            acceptance_criteria: 판정 기준 (예: {"agg_threshold": 5.0, "criteria_name": "Standard"})
        """
        try:
            # 1. 성공 여부 판정 (acceptance_criteria 기반)
            is_success = self._evaluate_success(assay_type, raw_data, acceptance_criteria)
            
            # 2. 신뢰도 점수 계산
            confidence_score = self._calculate_confidence(assay_type, raw_data, acceptance_criteria)
            
            # 3. 결과 저장
            result_id = str(uuid4())
            new_result = {
                "id": result_id,
                "run_id": run_id,
                "molecule_id": molecule_id,
                "assay_type": assay_type,
                "raw_data": raw_data,
                "is_success": is_success,
                "acceptance_criteria": acceptance_criteria,
                "confidence_score": confidence_score,
                "created_at": datetime.utcnow().isoformat()
            }
            
            res = supabase.table("assay_results").insert(new_result).execute()
            
            if res.data:
                logger.info(f"✅ Saved assay result {result_id}: {assay_type} -> {'Success' if is_success else 'Fail'}")
                return res.data[0]
            else:
                raise Exception("Failed to insert assay_result")
                
        except Exception as e:
            logger.error(f"Save assay result error: {e}")
            raise

    def _evaluate_success(
        self, 
        assay_type: str, 
        raw_data: Dict[str, Any], 
        criteria: Dict[str, Any]
    ) -> bool:
        """
        성공 여부 판정 (assay_type별 로직)
        """
        try:
            if assay_type == "aggregation":
                # Aggregation: raw_data의 값이 임계값보다 낮으면 성공
                agg_value = raw_data.get("aggregation_percent", 100)
                threshold = criteria.get("agg_threshold", 5.0)
                return agg_value < threshold
            
            elif assay_type == "binding":
                # Binding Affinity: raw_data의 값이 임계값보다 낮으면 성공 (낮을수록 좋음)
                kd_value = raw_data.get("kd_nm", 999)
                threshold = criteria.get("kd_threshold", 10.0)
                return kd_value <= threshold
            
            elif assay_type == "cytotoxicity":
                # IC50: 낮을수록 좋음
                ic50_value = raw_data.get("ic50_nm", 999)
                threshold = criteria.get("ic50_threshold", 100)
                return ic50_value <= threshold
            
            else:
                # 기본: raw_data에 "success" 키가 있으면 그 값 사용
                return raw_data.get("success", False)
                
        except Exception as e:
            logger.warning(f"Success evaluation error: {e}")
            return False

    def _calculate_confidence(
        self, 
        assay_type: str, 
        raw_data: Dict[str, Any], 
        criteria: Dict[str, Any]
    ) -> float:
        """
        Rule Confidence 점수 계산 (0.0 ~ 1.0)
        - 엄격한 기준에서 통과할수록 높은 신뢰도
        """
        try:
            base_score = 0.5  # 기본 점수
            
            if assay_type == "aggregation":
                agg_value = raw_data.get("aggregation_percent", 100)
                threshold = criteria.get("agg_threshold", 5.0)
                
                if agg_value < threshold:
                    # 임계값 대비 마진이 클수록 높은 점수
                    margin_ratio = (threshold - agg_value) / threshold
                    base_score = 0.5 + (margin_ratio * 0.5)
                else:
                    # 실패 시 낮은 점수
                    base_score = 0.3
            
            elif assay_type == "binding":
                kd_value = raw_data.get("kd_nm", 999)
                threshold = criteria.get("kd_threshold", 10.0)
                
                if kd_value <= threshold:
                    margin_ratio = (threshold - kd_value) / threshold
                    base_score = 0.5 + (margin_ratio * 0.5)
                else:
                    base_score = 0.3
            
            # 기준의 엄격도에 따른 보정 (엄격한 기준일수록 가중치 부여)
            strictness = criteria.get("strictness", 1.0)  # 1.0 = 표준, >1.0 = 엄격
            final_score = min(1.0, base_score * strictness)
            
            return round(final_score, 3)
            
        except Exception as e:
            logger.warning(f"Confidence calculation error: {e}")
            return 0.5

    async def get_results_by_run(self, run_id: str) -> List[Dict[str, Any]]:
        """특정 Design Run의 모든 Assay 결과 조회"""
        try:
            res = supabase.table("assay_results")\
                .select("*")\
                .eq("run_id", run_id)\
                .order("created_at", desc=True)\
                .execute()
            return res.data or []
        except Exception as e:
            logger.error(f"Get results by run error: {e}")
            return []

    async def get_results_by_molecule(self, molecule_id: str) -> List[Dict[str, Any]]:
        """특정 분자의 모든 Assay 결과 조회"""
        try:
            res = supabase.table("assay_results")\
                .select("*")\
                .eq("molecule_id", molecule_id)\
                .order("created_at", desc=True)\
                .execute()
            return res.data or []
        except Exception as e:
            logger.error(f"Get results by molecule error: {e}")
            return []


# Singleton
assay_result_service = AssayResultService()
