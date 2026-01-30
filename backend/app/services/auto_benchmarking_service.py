"""
Auto-Benchmarking Service
Catalog ADC와의 자동 비교 분석

Phase 2 Enhancement:
- catalog_adcs 테이블 연동
- Tanimoto 유사도 기반 매칭
- 비용/성능 비교 분석
"""
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field
from datetime import datetime
import logging

# RDKit imports (optional)
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, DataStructs
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

from app.core.supabase import get_supabase_client
from app.core.config import settings

logger = logging.getLogger(__name__)


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class CatalogADC:
    """카탈로그 ADC 정보"""
    id: str
    manufacturer: str
    catalog_number: str
    product_name: Optional[str]
    target: Optional[str]
    target_normalized: Optional[str]
    payload_smiles: Optional[str]
    linker_smiles: Optional[str]
    linker_type: Optional[str]
    dar: Optional[float]
    price_usd: Optional[float]
    availability: Optional[str]
    reported_kd_nm: Optional[float]
    half_life_hours: Optional[float]
    known_toxicity_score: Optional[float]


@dataclass
class SimilarityResult:
    """유사도 분석 결과"""
    structural_similarity: float  # Tanimoto (0-1)
    target_match: bool
    linker_type_match: bool
    overall_score: float  # 가중 평균 (0-100)


@dataclass
class ComparisonMetric:
    """비교 지표"""
    name: str
    our_value: Optional[float]
    catalog_value: Optional[float]
    improvement_pct: float
    is_favorable: bool
    confidence: float


@dataclass
class BenchmarkComparison:
    """벤치마크 비교 결과"""
    catalog_adc: CatalogADC
    similarity: SimilarityResult
    comparison_metrics: Dict[str, ComparisonMetric]
    verdict: Dict[str, Any]
    reasoning: str


# ============================================================================
# Auto-Benchmarking Service
# ============================================================================

class AutoBenchmarkingService:
    """
    Catalog ADC 자동 벤치마킹 서비스

    신규 설계안과 기존 상용 ADC 시약의 자동 비교 분석
    """

    # 가중치 설정
    SIMILARITY_WEIGHTS = {
        "structural": 0.4,  # Tanimoto 유사도
        "target": 0.3,      # 타겟 일치
        "linker": 0.2,      # 링커 타입 일치
        "dar": 0.1          # DAR 근접성
    }

    def __init__(self):
        self.supabase = get_supabase_client()

    async def auto_benchmark(
        self,
        target_antigen: str,
        payload_smiles: str,
        linker_type: Optional[str] = None,
        linker_smiles: Optional[str] = None,
        our_metrics: Optional[Dict[str, Any]] = None,
        dar: Optional[float] = None,
        limit: int = 10
    ) -> Dict[str, Any]:
        """
        자동 벤치마킹 수행

        Args:
            target_antigen: 타겟 항원
            payload_smiles: 페이로드 SMILES
            linker_type: 링커 타입 (cleavable, non-cleavable)
            linker_smiles: 링커 SMILES
            our_metrics: 우리 ADC 메트릭
            dar: DAR 값
            limit: 검색 제한

        Returns:
            벤치마크 비교 결과
        """
        our_metrics = our_metrics or {}

        # 1. 유사 Catalog ADC 검색
        similar_products = await self._search_catalog_adcs(
            target=target_antigen,
            linker_type=linker_type,
            limit=limit
        )

        if not similar_products:
            return {
                "success": False,
                "message": "No similar catalog ADCs found",
                "benchmarks": []
            }

        # 2. 구조적 유사도 계산 및 점수화
        scored_products = []
        our_fp = self._get_fingerprint(payload_smiles) if payload_smiles else None

        for product in similar_products:
            catalog_adc = self._dict_to_catalog_adc(product)

            # 유사도 계산
            similarity = self._calculate_similarity(
                our_fp=our_fp,
                catalog_smiles=catalog_adc.payload_smiles or catalog_adc.linker_smiles,
                target_match=(
                    target_antigen and catalog_adc.target_normalized and
                    target_antigen.upper() in catalog_adc.target_normalized.upper()
                ),
                linker_match=(
                    linker_type and catalog_adc.linker_type and
                    linker_type.lower() == catalog_adc.linker_type.lower()
                ),
                our_dar=dar,
                catalog_dar=catalog_adc.dar
            )

            scored_products.append({
                "catalog_adc": catalog_adc,
                "similarity": similarity
            })

        # 3. 상위 3개 선택
        sorted_products = sorted(
            scored_products,
            key=lambda x: x["similarity"].overall_score,
            reverse=True
        )[:3]

        # 4. 상세 비교 분석
        benchmarks = []
        for item in sorted_products:
            comparison = await self._compare_with_catalog(
                our_metrics=our_metrics,
                catalog_adc=item["catalog_adc"],
                similarity=item["similarity"]
            )
            benchmarks.append(comparison)

        # 5. 종합 분석
        overall_analysis = self._generate_overall_analysis(benchmarks, target_antigen)

        return {
            "success": True,
            "target": target_antigen,
            "total_similar_found": len(similar_products),
            "benchmarks": [self._benchmark_to_dict(b) for b in benchmarks],
            "overall_analysis": overall_analysis,
            "recommendation": self._get_overall_recommendation(benchmarks)
        }

    async def _search_catalog_adcs(
        self,
        target: Optional[str] = None,
        linker_type: Optional[str] = None,
        limit: int = 10
    ) -> List[Dict[str, Any]]:
        """
        Catalog ADC 검색

        Args:
            target: 타겟 항원
            linker_type: 링커 타입
            limit: 검색 제한

        Returns:
            매칭 카탈로그 ADC 목록
        """
        try:
            query = self.supabase.table("catalog_adcs").select(
                "id, manufacturer, catalog_number, product_name, "
                "target, target_normalized, payload_smiles, linker_smiles, "
                "linker_type, dar, price_usd, availability, "
                "reported_kd_nm, half_life_hours, known_toxicity_score"
            )

            # 타겟 기반 필터링
            if target:
                query = query.or_(
                    f"target_normalized.ilike.%{target}%,target.ilike.%{target}%"
                )

            # 링커 타입 필터링
            if linker_type:
                query = query.eq("linker_type", linker_type)

            # 가격이 있는 것 우선
            result = query.not_.is_("price_usd", "null").limit(limit).execute()

            products = result.data or []

            # 추가 검색 (가격 없는 것도 포함)
            if len(products) < limit:
                remaining = limit - len(products)
                existing_ids = [p["id"] for p in products]

                extra_query = self.supabase.table("catalog_adcs").select(
                    "id, manufacturer, catalog_number, product_name, "
                    "target, target_normalized, payload_smiles, linker_smiles, "
                    "linker_type, dar, price_usd, availability, "
                    "reported_kd_nm, half_life_hours, known_toxicity_score"
                )

                if target:
                    extra_query = extra_query.or_(
                        f"target_normalized.ilike.%{target}%,target.ilike.%{target}%"
                    )

                extra_result = extra_query.limit(remaining + len(existing_ids)).execute()

                for p in (extra_result.data or []):
                    if p["id"] not in existing_ids:
                        products.append(p)
                        if len(products) >= limit:
                            break

            logger.info(f"[AutoBenchmarking] Found {len(products)} catalog ADCs for target: {target}")
            return products

        except Exception as e:
            logger.error(f"[AutoBenchmarking] Search error: {e}")
            return []

    def _get_fingerprint(self, smiles: str) -> Optional[Any]:
        """SMILES에서 Morgan Fingerprint 생성"""
        if not RDKIT_AVAILABLE or not smiles:
            return None

        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            return None
        except Exception as e:
            logger.warning(f"[AutoBenchmarking] Fingerprint error: {e}")
            return None

    def _calculate_tanimoto(self, fp1: Any, fp2: Any) -> float:
        """Tanimoto 유사도 계산"""
        if not RDKIT_AVAILABLE or fp1 is None or fp2 is None:
            return 0.5  # Default similarity

        try:
            return DataStructs.TanimotoSimilarity(fp1, fp2)
        except Exception:
            return 0.5

    def _calculate_similarity(
        self,
        our_fp: Optional[Any],
        catalog_smiles: Optional[str],
        target_match: bool,
        linker_match: bool,
        our_dar: Optional[float],
        catalog_dar: Optional[float]
    ) -> SimilarityResult:
        """
        종합 유사도 계산

        Returns:
            SimilarityResult with overall_score (0-100)
        """
        # 구조적 유사도
        catalog_fp = self._get_fingerprint(catalog_smiles) if catalog_smiles else None
        structural_sim = self._calculate_tanimoto(our_fp, catalog_fp)

        # DAR 유사도 (근접성)
        dar_sim = 1.0
        if our_dar is not None and catalog_dar is not None:
            dar_diff = abs(our_dar - catalog_dar)
            dar_sim = max(0, 1 - (dar_diff / 4))  # DAR 4 차이면 0

        # 가중 점수 계산
        overall = (
            self.SIMILARITY_WEIGHTS["structural"] * structural_sim * 100 +
            self.SIMILARITY_WEIGHTS["target"] * (100 if target_match else 30) +
            self.SIMILARITY_WEIGHTS["linker"] * (100 if linker_match else 50) +
            self.SIMILARITY_WEIGHTS["dar"] * dar_sim * 100
        )

        return SimilarityResult(
            structural_similarity=round(structural_sim, 3),
            target_match=target_match,
            linker_type_match=linker_match,
            overall_score=round(overall, 1)
        )

    async def _compare_with_catalog(
        self,
        our_metrics: Dict[str, Any],
        catalog_adc: CatalogADC,
        similarity: SimilarityResult
    ) -> BenchmarkComparison:
        """
        개별 카탈로그 ADC와 상세 비교

        Args:
            our_metrics: 우리 ADC 메트릭
            catalog_adc: 카탈로그 ADC 정보
            similarity: 유사도 결과

        Returns:
            BenchmarkComparison 결과
        """
        comparison_metrics = {}

        # 1. 독성 비교
        our_tox = our_metrics.get("predicted_toxicity", 50)
        catalog_tox = catalog_adc.known_toxicity_score or 60
        tox_improvement = ((catalog_tox - our_tox) / catalog_tox) * 100 if catalog_tox else 0
        comparison_metrics["toxicity"] = ComparisonMetric(
            name="Predicted Toxicity",
            our_value=our_tox,
            catalog_value=catalog_tox,
            improvement_pct=round(tox_improvement, 1),
            is_favorable=tox_improvement > 0,
            confidence=0.7
        )

        # 2. 결합 친화도 비교 (Kd)
        our_kd = our_metrics.get("predicted_kd_nm", 10)
        catalog_kd = catalog_adc.reported_kd_nm or 15
        kd_improvement = (catalog_kd / our_kd) if our_kd > 0 else 1.0
        comparison_metrics["binding_affinity"] = ComparisonMetric(
            name="Binding Affinity (Kd)",
            our_value=our_kd,
            catalog_value=catalog_kd,
            improvement_pct=round((kd_improvement - 1) * 100, 1),
            is_favorable=kd_improvement > 1,
            confidence=0.6
        )

        # 3. 안정성 비교 (반감기)
        our_half_life = our_metrics.get("half_life_hours", 72)
        catalog_half_life = catalog_adc.half_life_hours or 48
        stability_improvement = ((our_half_life - catalog_half_life) / catalog_half_life) * 100 if catalog_half_life else 0
        comparison_metrics["stability"] = ComparisonMetric(
            name="Half-life (hours)",
            our_value=our_half_life,
            catalog_value=catalog_half_life,
            improvement_pct=round(stability_improvement, 1),
            is_favorable=stability_improvement > 0,
            confidence=0.65
        )

        # 4. 합성 용이성 비교
        our_sa = our_metrics.get("sa_score", 4.0)
        catalog_sa = 4.5  # 기본 추정값
        comparison_metrics["synthesizability"] = ComparisonMetric(
            name="Synthetic Accessibility",
            our_value=our_sa,
            catalog_value=catalog_sa,
            improvement_pct=round(((catalog_sa - our_sa) / catalog_sa) * 100, 1),
            is_favorable=our_sa < catalog_sa,
            confidence=0.8
        )

        # 5. 비용 비교
        catalog_price = catalog_adc.price_usd or 3000
        estimated_synthesis = self._estimate_synthesis_cost(our_metrics)
        cost_saving = ((catalog_price - estimated_synthesis) / catalog_price) * 100 if catalog_price else 0
        comparison_metrics["cost"] = ComparisonMetric(
            name="Cost (USD)",
            our_value=estimated_synthesis,
            catalog_value=catalog_price,
            improvement_pct=round(cost_saving, 1),
            is_favorable=cost_saving > 0,
            confidence=0.5  # 비용 추정은 불확실
        )

        # 6. 종합 판정
        verdict = self._generate_verdict(comparison_metrics, similarity)

        # 7. 추론 생성
        reasoning = self._generate_reasoning(comparison_metrics, catalog_adc, similarity)

        return BenchmarkComparison(
            catalog_adc=catalog_adc,
            similarity=similarity,
            comparison_metrics=comparison_metrics,
            verdict=verdict,
            reasoning=reasoning
        )

    def _estimate_synthesis_cost(self, our_metrics: Dict[str, Any]) -> float:
        """합성 비용 추정"""
        base_cost = 1000  # 기본 비용

        sa_score = our_metrics.get("sa_score", 4.0)
        # SA 점수가 높을수록 비용 증가
        sa_multiplier = 1 + (sa_score - 3) * 0.3

        mw = our_metrics.get("mw", 500)
        # 분자량이 클수록 비용 증가
        mw_multiplier = 1 + (mw - 400) / 1000

        return round(base_cost * sa_multiplier * mw_multiplier, 2)

    def _generate_verdict(
        self,
        comparison_metrics: Dict[str, ComparisonMetric],
        similarity: SimilarityResult
    ) -> Dict[str, Any]:
        """종합 판정 생성"""
        advantages = []
        caveats = []

        for key, metric in comparison_metrics.items():
            if metric.is_favorable and metric.improvement_pct > 10:
                if key == "toxicity":
                    advantages.append(f"{metric.improvement_pct:.0f}% lower predicted toxicity")
                elif key == "binding_affinity":
                    advantages.append(f"{metric.improvement_pct:.0f}% better binding affinity")
                elif key == "stability":
                    advantages.append(f"{metric.improvement_pct:.0f}% longer half-life")
                elif key == "cost":
                    advantages.append(f"{metric.improvement_pct:.0f}% cost reduction potential")

            if not metric.is_favorable and abs(metric.improvement_pct) > 20:
                caveats.append(f"{metric.name}: {abs(metric.improvement_pct):.0f}% lower than catalog")

        # 추천 결정
        if len(advantages) >= 2:
            recommendation = "use_our_design"
            primary_advantage = advantages[0] if advantages else "Novel design"
        elif len(advantages) == 1:
            recommendation = "use_our_design"
            primary_advantage = advantages[0]
        elif similarity.overall_score > 80:
            recommendation = "similar_performance"
            primary_advantage = "Comparable to existing products"
        else:
            recommendation = "consider_catalog"
            primary_advantage = "Catalog product may be more cost-effective"

        return {
            "recommendation": recommendation,
            "primary_advantage": primary_advantage,
            "advantages": advantages,
            "caveats": caveats
        }

    def _generate_reasoning(
        self,
        comparison_metrics: Dict[str, ComparisonMetric],
        catalog_adc: CatalogADC,
        similarity: SimilarityResult
    ) -> str:
        """비교 논리 생성"""
        parts = [
            f"[vs {catalog_adc.manufacturer} {catalog_adc.catalog_number}]"
        ]

        parts.append(f"구조적 유사도: {similarity.structural_similarity*100:.1f}% (Tanimoto)")
        parts.append(f"타겟 일치: {'Yes' if similarity.target_match else 'No'}")
        parts.append(f"링커 타입 일치: {'Yes' if similarity.linker_type_match else 'No'}")
        parts.append(f"종합 유사도 점수: {similarity.overall_score:.1f}/100")

        for key, metric in comparison_metrics.items():
            direction = "우수" if metric.is_favorable else "열위"
            parts.append(
                f"• {metric.name}: {metric.our_value} vs {metric.catalog_value} "
                f"({direction}, {metric.improvement_pct:+.1f}%)"
            )

        if catalog_adc.price_usd:
            parts.append(f"카탈로그 가격: ${catalog_adc.price_usd:,.2f}")
        if catalog_adc.availability:
            parts.append(f"구매 가능성: {catalog_adc.availability}")

        return "\n".join(parts)

    def _generate_overall_analysis(
        self,
        benchmarks: List[BenchmarkComparison],
        target: str
    ) -> Dict[str, Any]:
        """전체 분석 요약"""
        if not benchmarks:
            return {
                "summary": "No catalog ADCs available for comparison",
                "market_position": "novel"
            }

        # 평균 유사도
        avg_similarity = sum(b.similarity.overall_score for b in benchmarks) / len(benchmarks)

        # 우수 지표 카운트
        advantage_counts = {}
        for b in benchmarks:
            for key, metric in b.comparison_metrics.items():
                if metric.is_favorable:
                    advantage_counts[key] = advantage_counts.get(key, 0) + 1

        # 가장 빈번한 우위 지표
        top_advantages = sorted(advantage_counts.items(), key=lambda x: x[1], reverse=True)[:3]

        # 시장 포지션
        if avg_similarity > 80:
            market_position = "me_too"
            position_desc = "기존 제품과 유사한 프로파일"
        elif avg_similarity > 50:
            market_position = "differentiated"
            position_desc = "기존 제품 대비 차별화된 특성"
        else:
            market_position = "novel"
            position_desc = "시장에 유사 제품이 적은 혁신적 설계"

        return {
            "summary": f"{target} 타겟 Catalog ADC {len(benchmarks)}종과 비교 분석 완료",
            "average_similarity": round(avg_similarity, 1),
            "market_position": market_position,
            "position_description": position_desc,
            "top_advantages": [
                {"metric": k, "count": v, "percentage": round(v/len(benchmarks)*100, 1)}
                for k, v in top_advantages
            ],
            "total_comparisons": len(benchmarks)
        }

    def _get_overall_recommendation(
        self,
        benchmarks: List[BenchmarkComparison]
    ) -> Dict[str, Any]:
        """전체 추천"""
        if not benchmarks:
            return {
                "action": "proceed",
                "confidence": 0.5,
                "reason": "No comparable catalog products found"
            }

        # 추천 집계
        use_design = sum(1 for b in benchmarks if b.verdict["recommendation"] == "use_our_design")
        consider_catalog = sum(1 for b in benchmarks if b.verdict["recommendation"] == "consider_catalog")

        if use_design >= len(benchmarks) * 0.6:
            action = "proceed_with_confidence"
            confidence = 0.85
            reason = f"{use_design}/{len(benchmarks)} 비교에서 우리 설계가 우위"
        elif use_design >= len(benchmarks) * 0.3:
            action = "proceed_with_optimization"
            confidence = 0.65
            reason = "일부 지표에서 개선 필요하나 진행 권고"
        else:
            action = "review_alternatives"
            confidence = 0.4
            reason = "기존 카탈로그 제품 검토 권고"

        return {
            "action": action,
            "confidence": confidence,
            "reason": reason,
            "use_design_count": use_design,
            "consider_catalog_count": consider_catalog,
            "total_comparisons": len(benchmarks)
        }

    # =========================================================================
    # Helper Methods
    # =========================================================================

    def _dict_to_catalog_adc(self, data: Dict[str, Any]) -> CatalogADC:
        """Dict를 CatalogADC로 변환"""
        return CatalogADC(
            id=data.get("id", ""),
            manufacturer=data.get("manufacturer", "Unknown"),
            catalog_number=data.get("catalog_number", ""),
            product_name=data.get("product_name"),
            target=data.get("target"),
            target_normalized=data.get("target_normalized"),
            payload_smiles=data.get("payload_smiles"),
            linker_smiles=data.get("linker_smiles"),
            linker_type=data.get("linker_type"),
            dar=data.get("dar"),
            price_usd=data.get("price_usd"),
            availability=data.get("availability"),
            reported_kd_nm=data.get("reported_kd_nm"),
            half_life_hours=data.get("half_life_hours"),
            known_toxicity_score=data.get("known_toxicity_score")
        )

    def _benchmark_to_dict(self, benchmark: BenchmarkComparison) -> Dict[str, Any]:
        """BenchmarkComparison을 Dict로 변환"""
        return {
            "catalog_adc": {
                "id": benchmark.catalog_adc.id,
                "manufacturer": benchmark.catalog_adc.manufacturer,
                "catalog_number": benchmark.catalog_adc.catalog_number,
                "product_name": benchmark.catalog_adc.product_name,
                "target": benchmark.catalog_adc.target,
                "linker_type": benchmark.catalog_adc.linker_type,
                "dar": benchmark.catalog_adc.dar,
                "price_usd": benchmark.catalog_adc.price_usd,
                "availability": benchmark.catalog_adc.availability
            },
            "similarity": {
                "structural_similarity": benchmark.similarity.structural_similarity,
                "target_match": benchmark.similarity.target_match,
                "linker_type_match": benchmark.similarity.linker_type_match,
                "overall_score": benchmark.similarity.overall_score
            },
            "comparison_metrics": {
                key: {
                    "name": metric.name,
                    "our_value": metric.our_value,
                    "catalog_value": metric.catalog_value,
                    "improvement_pct": metric.improvement_pct,
                    "is_favorable": metric.is_favorable,
                    "confidence": metric.confidence
                }
                for key, metric in benchmark.comparison_metrics.items()
            },
            "verdict": benchmark.verdict,
            "reasoning": benchmark.reasoning
        }


# ============================================================================
# Module-level Functions
# ============================================================================

_service_instance = None


def get_auto_benchmarking_service() -> AutoBenchmarkingService:
    """싱글톤 인스턴스 반환"""
    global _service_instance
    if _service_instance is None:
        _service_instance = AutoBenchmarkingService()
    return _service_instance


async def run_auto_benchmark(
    target_antigen: str,
    payload_smiles: str,
    linker_type: Optional[str] = None,
    our_metrics: Optional[Dict[str, Any]] = None,
    dar: Optional[float] = None
) -> Dict[str, Any]:
    """자동 벤치마킹 실행 (편의 함수)"""
    service = get_auto_benchmarking_service()
    return await service.auto_benchmark(
        target_antigen=target_antigen,
        payload_smiles=payload_smiles,
        linker_type=linker_type,
        our_metrics=our_metrics,
        dar=dar
    )
