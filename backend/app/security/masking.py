"""
Data Masking Service - Backend Level Masking
Constraint 3: Free Tier 데이터 백엔드 레벨 마스킹
"""
from typing import Any, Dict, List
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class MaskingLevel(Enum):
    """마스킹 레벨"""
    NONE = "none"           # Premium: 전체 공개
    PARTIAL = "partial"     # Free: 부분 마스킹
    FULL = "full"           # 완전 마스킹


class DataMaskingService:
    """
    Constraint 3 준수: 백엔드 API 레벨에서 데이터 마스킹

    프론트엔드가 아닌 서버에서 직접 마스킹하여 보안 강화
    Free Tier 사용자에게 핵심 데이터를 숨기고 Premium 업그레이드 유도
    """

    def mask_session_result(self, data: Dict, tier: str) -> Dict:
        """
        세션 결과 마스킹

        Args:
            data: 원본 세션 데이터
            tier: 사용자 티어 ('free' or 'premium')

        Returns:
            마스킹된 데이터
        """
        if tier == "premium":
            return data  # Premium: 전체 공개

        # Free tier: 핵심 데이터 마스킹
        masked = data.copy()

        # 1. SMILES 마스킹 (첫 10자만 표시)
        if "current_smiles" in masked and masked["current_smiles"]:
            smiles = masked["current_smiles"]
            masked["current_smiles"] = self.mask_smiles(smiles, MaskingLevel.PARTIAL)
            masked["smiles_masked"] = True

        # 2. 후보 리스트 마스킹 (1개만 표시, 나머지 블러)
        if "candidates" in masked and masked["candidates"]:
            candidates = masked["candidates"]
            if len(candidates) > 1:
                masked["candidates"] = [
                    candidates[0],  # 첫 번째만 표시
                    *[self._mask_candidate(c) for c in candidates[1:]]
                ]
                masked["candidates_limited"] = True
                masked["total_candidates_available"] = len(candidates)

        # 3. 점수 마스킹
        if "calculated_metrics" in masked and masked["calculated_metrics"]:
            masked["calculated_metrics"] = self._mask_metrics(
                masked["calculated_metrics"], tier
            )

        # 4. 리포트 마스킹
        if "final_report" in masked and masked["final_report"]:
            masked["final_report"] = self._mask_report(masked["final_report"])

        return masked

    def _mask_candidate(self, candidate: Dict) -> Dict:
        """개별 후보 마스킹"""
        return {
            "rank": candidate.get("rank"),
            "smiles": "********** (Premium Only)",
            "score": "?.??",
            "metrics": None,
            "validation": None,
            "is_masked": True
        }

    def _mask_metrics(self, metrics: Dict, tier: str) -> Dict:
        """메트릭 마스킹"""
        if tier == "premium":
            return metrics

        return {
            "mw": metrics.get("mw"),  # MW는 공개
            "logp": "***" if metrics.get("logp") else None,
            "hbd": metrics.get("hbd"),  # 기본 정보는 공개
            "hba": metrics.get("hba"),
            "sa_score": "Premium Only",
            "tpsa": "Premium Only",
            "is_masked": True
        }

    def _mask_report(self, report: Dict) -> Dict:
        """리포트 마스킹"""
        summary = report.get("summary", {})
        return {
            "summary": {
                "target": summary.get("target"),
                "indication": summary.get("indication"),
                "smiles_generated": summary.get("smiles_generated"),
                "candidates_count": summary.get("candidates_count")
            },
            "validation": {
                "chemistry": "Premium Only",
                "constraints": {
                    "passed": report.get("validation", {}).get("constraints", {}).get("passed"),
                    "details": "Premium Only"
                }
            },
            "risk_assessment": {
                "level": report.get("risk_assessment", {}).get("level"),
                "score": "?/10"
            },
            "decision": {
                "approved": report.get("decision", {}).get("approved"),
                "reasoning": report.get("decision", {}).get("reasoning", "")[:100] + "... (Premium Only)",
                "confidence": "Premium Only"
            },
            "pmid_references": (report.get("pmid_references") or [])[:1],  # 1개만 표시
            "is_masked": True
        }

    def mask_smiles(self, smiles: str, level: MaskingLevel) -> str:
        """
        SMILES 문자열 마스킹

        Args:
            smiles: 원본 SMILES
            level: 마스킹 레벨

        Returns:
            마스킹된 SMILES
        """
        if not smiles:
            return smiles

        if level == MaskingLevel.NONE:
            return smiles
        elif level == MaskingLevel.PARTIAL:
            # 첫 10자 + 마스크
            visible = smiles[:10]
            masked_length = min(20, len(smiles) - 10)
            return f"{visible}...{'*' * masked_length}"
        else:
            return "*" * min(30, len(smiles))

    def mask_list_for_free_tier(
        self,
        items: List[Dict],
        visible_count: int = 1,
        mask_fields: List[str] = None
    ) -> Dict:
        """
        리스트 데이터 마스킹 (Free Tier용)

        Args:
            items: 원본 리스트
            visible_count: 보여줄 아이템 수
            mask_fields: 마스킹할 필드들

        Returns:
            {
                "items": 마스킹된 리스트,
                "total_available": 전체 개수,
                "is_limited": True
            }
        """
        if not items:
            return {"items": [], "total_available": 0, "is_limited": False}

        mask_fields = mask_fields or ["smiles", "structure", "score"]

        visible_items = items[:visible_count]
        masked_items = []

        for item in items[visible_count:]:
            masked_item = item.copy()
            for field in mask_fields:
                if field in masked_item:
                    masked_item[field] = "Premium Only"
            masked_item["is_masked"] = True
            masked_items.append(masked_item)

        return {
            "items": visible_items + masked_items,
            "total_available": len(items),
            "is_limited": len(items) > visible_count
        }
