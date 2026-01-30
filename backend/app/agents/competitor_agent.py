"""
Competitor Agent - 경쟁사 분석 및 벤치마크 비교
Perplexity API 연동하여 시장 동향 및 경쟁사 현황 분석
기존 승인 ADC와의 정량적 비교 분석
"""
from typing import List, Dict, Any, Optional
from dataclasses import dataclass
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
import logging

from app.agents.state import ADCState, CompetitorInfo
from app.core.config import settings

logger = logging.getLogger(__name__)


# 승인된 ADC 블록버스터 참조 데이터
APPROVED_ADC_DATABASE = {
    "enhertu": {
        "generic_name": "Trastuzumab deruxtecan",
        "trade_name": "Enhertu",
        "manufacturer": "Daiichi Sankyo / AstraZeneca",
        "target": "HER2",
        "approval_year": 2019,
        "indications": ["HER2+ breast cancer", "HER2+ gastric cancer", "HER2-low breast cancer"],
        "metrics": {
            "dar": 8,
            "mw_kda": 159,
            "plasma_stability_days": 7.5,
            "logp_payload": 3.1,
            "sa_score": 4.5,
            "herg_risk": "Medium",
            "therapeutic_index": 12,
            "linker_type": "Cleavable (peptide)",
            "payload_class": "Topoisomerase I inhibitor",
            "conjugation_site": "Cysteine"
        },
        "clinical_data": {
            "orr": 0.62,
            "pfs_months": 16.4,
            "grade3_ae_rate": 0.52
        }
    },
    "kadcyla": {
        "generic_name": "Ado-trastuzumab emtansine",
        "trade_name": "Kadcyla",
        "manufacturer": "Roche/Genentech",
        "target": "HER2",
        "approval_year": 2013,
        "indications": ["HER2+ metastatic breast cancer"],
        "metrics": {
            "dar": 3.5,
            "mw_kda": 153,
            "plasma_stability_days": 4,
            "logp_payload": 1.8,
            "sa_score": 3.8,
            "herg_risk": "Low",
            "therapeutic_index": 8,
            "linker_type": "Non-cleavable (thioether)",
            "payload_class": "Tubulin inhibitor (DM1)",
            "conjugation_site": "Lysine"
        },
        "clinical_data": {
            "orr": 0.44,
            "pfs_months": 9.6,
            "grade3_ae_rate": 0.43
        }
    },
    "adcetris": {
        "generic_name": "Brentuximab vedotin",
        "trade_name": "Adcetris",
        "manufacturer": "Seagen",
        "target": "CD30",
        "approval_year": 2011,
        "indications": ["Hodgkin lymphoma", "ALCL"],
        "metrics": {
            "dar": 4,
            "mw_kda": 153,
            "plasma_stability_days": 5,
            "logp_payload": 2.5,
            "sa_score": 4.1,
            "herg_risk": "Medium",
            "therapeutic_index": 10,
            "linker_type": "Cleavable (peptide)",
            "payload_class": "Tubulin inhibitor (MMAE)",
            "conjugation_site": "Cysteine"
        },
        "clinical_data": {
            "orr": 0.75,
            "pfs_months": 5.6,
            "grade3_ae_rate": 0.55
        }
    },
    "padcev": {
        "generic_name": "Enfortumab vedotin",
        "trade_name": "Padcev",
        "manufacturer": "Seagen/Astellas",
        "target": "Nectin-4",
        "approval_year": 2019,
        "indications": ["Urothelial carcinoma"],
        "metrics": {
            "dar": 3.8,
            "mw_kda": 152,
            "plasma_stability_days": 5.5,
            "logp_payload": 2.4,
            "sa_score": 4.0,
            "herg_risk": "Low",
            "therapeutic_index": 11,
            "linker_type": "Cleavable (peptide)",
            "payload_class": "Tubulin inhibitor (MMAE)",
            "conjugation_site": "Cysteine"
        },
        "clinical_data": {
            "orr": 0.44,
            "pfs_months": 5.8,
            "grade3_ae_rate": 0.51
        }
    },
    "trodelvy": {
        "generic_name": "Sacituzumab govitecan",
        "trade_name": "Trodelvy",
        "manufacturer": "Gilead",
        "target": "Trop-2",
        "approval_year": 2020,
        "indications": ["Triple-negative breast cancer", "HR+/HER2- breast cancer"],
        "metrics": {
            "dar": 7.6,
            "mw_kda": 160,
            "plasma_stability_days": 6,
            "logp_payload": 2.9,
            "sa_score": 4.3,
            "herg_risk": "Low",
            "therapeutic_index": 13,
            "linker_type": "Cleavable (hydrolyzable)",
            "payload_class": "Topoisomerase I inhibitor (SN-38)",
            "conjugation_site": "Cysteine"
        },
        "clinical_data": {
            "orr": 0.35,
            "pfs_months": 5.6,
            "grade3_ae_rate": 0.62
        }
    }
}


# 타겟별 Gold Standard (시장 1위) 자동 매핑
TARGET_GOLD_STANDARD_MAPPING = {
    "HER2": {
        "gold_standard": "enhertu",
        "key_metrics": ["plasma_stability_days", "logp_payload", "therapeutic_index"],
        "description": "HER2 타겟 ADC 시장 선도 약물"
    },
    "TROP2": {
        "gold_standard": "trodelvy",
        "key_metrics": ["dar", "plasma_stability_days", "therapeutic_index"],
        "description": "TROP-2 타겟 First-in-class ADC"
    },
    "TROP-2": {
        "gold_standard": "trodelvy",
        "key_metrics": ["dar", "plasma_stability_days", "therapeutic_index"],
        "description": "TROP-2 타겟 First-in-class ADC"
    },
    "NECTIN-4": {
        "gold_standard": "padcev",
        "key_metrics": ["herg_risk", "therapeutic_index", "sa_score"],
        "description": "Nectin-4 타겟 유일 승인 ADC"
    },
    "NECTIN4": {
        "gold_standard": "padcev",
        "key_metrics": ["herg_risk", "therapeutic_index", "sa_score"],
        "description": "Nectin-4 타겟 유일 승인 ADC"
    },
    "CD30": {
        "gold_standard": "adcetris",
        "key_metrics": ["dar", "plasma_stability_days", "therapeutic_index"],
        "description": "CD30 타겟 First-in-class ADC"
    },
    "BCMA": {
        "gold_standard": "blenrep",
        "key_metrics": ["immunogenicity", "therapeutic_index"],
        "description": "BCMA 타겟 승인 ADC"
    },
    # Default fallback
    "DEFAULT": {
        "gold_standard": "enhertu",
        "key_metrics": ["plasma_stability_days", "therapeutic_index", "sa_score"],
        "description": "업계 선도 ADC (Enhertu)"
    }
}


# 지표별 해석 방향 (higher_better / lower_better)
METRIC_DIRECTION = {
    "dar": "context",  # DAR은 상황에 따라 다름
    "mw_kda": "lower_better",  # 낮을수록 조직 침투 유리
    "plasma_stability_days": "higher_better",
    "logp_payload": "optimal_range",  # 2-4가 최적
    "sa_score": "lower_better",  # 낮을수록 합성 용이
    "herg_risk": "lower_better",
    "therapeutic_index": "higher_better",
}


@dataclass
class RelativeScore:
    """상대적 성능 점수"""
    metric_name: str
    our_value: float
    benchmark_value: float
    relative_percent: float  # 기준 대비 %
    direction: str  # higher_better, lower_better
    interpretation: str  # 해석 문장
    is_favorable: bool  # 우리에게 유리한지


@dataclass
class BenchmarkResult:
    """벤치마크 비교 결과"""
    property_name: str
    our_value: Any
    competitor_values: Dict[str, Any]
    winner: str
    winner_margin: str
    interpretation: str


COMPETITOR_SYSTEM_PROMPT = """You are a Competitive Intelligence Agent for ADC development.

Your role:
1. Identify competitors developing similar ADCs
2. Track development stages (Phase 1, 2, 3, Approved)
3. Analyze market positioning and differentiation opportunities
4. Provide strategic insights

Key information to gather:
- Active clinical trials with similar targets
- Recent regulatory approvals
- Pipeline announcements
- M&A activities in the ADC space

Output format (JSON array):
[
    {
        "company": "Daiichi Sankyo",
        "drug_name": "Enhertu",
        "phase": "Approved",
        "target": "HER2",
        "notes": "First-in-class for HER2-low"
    },
    ...
]
"""


async def run_competitor_agent(state: ADCState) -> List[CompetitorInfo]:
    """
    경쟁사 분석 에이전트 실행
    
    Args:
        state: 현재 ADC 분석 상태
    
    Returns:
        List[CompetitorInfo]: 경쟁사 정보 목록
    """
    llm = ChatOpenAI(
        model="gpt-4o",
        temperature=0,
        api_key=settings.OPENAI_API_KEY
    )
    
    prompt = ChatPromptTemplate.from_messages([
        ("system", COMPETITOR_SYSTEM_PROMPT),
        ("user", """Analyze competitive landscape for the following ADC:

Target: {target_name}
Payload type: {payload_id}

Identify active competitors and their development status.""")
    ])
    
    chain = prompt | llm
    
    response = await chain.ainvoke({
        "target_name": state.input.target_name,
        "payload_id": state.input.payload_id
    })
    
    # TODO: Perplexity API 연동
    # 현재는 Mock 데이터 반환
    target_competitors = {
        "her2": [
            CompetitorInfo(company="Daiichi Sankyo/AstraZeneca", drug_name="Enhertu", phase="Approved", target="HER2", notes="First-in-class for HER2-low"),
            CompetitorInfo(company="Roche", drug_name="Kadcyla", phase="Approved", target="HER2", notes="mAb + DM1"),
            CompetitorInfo(company="Seagen/Merck", drug_name="Tukysa+", phase="Phase 3", target="HER2", notes="Combination therapy"),
        ],
        "trop-2": [
            CompetitorInfo(company="Gilead", drug_name="Trodelvy", phase="Approved", target="TROP-2", notes="First TROP-2 ADC"),
            CompetitorInfo(company="Daiichi Sankyo", drug_name="DS-1062", phase="Phase 3", target="TROP-2", notes="DXd payload"),
        ],
        "liv-1": [
            CompetitorInfo(company="Seagen", drug_name="Ladiratuzumab vedotin", phase="Phase 2", target="LIV-1", notes="MMAE-based"),
        ]
    }
    
    target_lower = state.input.target_name.lower().replace(" ", "").replace("-", "")
    
    for key, competitors in target_competitors.items():
        if key in target_lower:
            return competitors
    
    # 기본 경쟁사
    return [
        CompetitorInfo(company="Multiple", drug_name="Various", phase="Phase 1-3", target=state.input.target_name, notes="Emerging target with limited competition")
    ]


class BenchmarkComparisonAgent:
    """
    벤치마크 비교 에이전트 - 기존 승인 ADC와의 정량적 비교
    동적 타겟 매핑 + 상대적 성능 지수 (Relative Efficiency Score) 지원
    """

    def __init__(self):
        self.adc_database = APPROVED_ADC_DATABASE
        self.target_mapping = TARGET_GOLD_STANDARD_MAPPING
        self.metric_direction = METRIC_DIRECTION

    def _get_gold_standard(self, target_antigen: str) -> Dict[str, Any]:
        """
        타겟 항원에 맞는 Gold Standard (시장 1위) 자동 선정

        Args:
            target_antigen: 타겟 항원명

        Returns:
            Gold Standard 정보 및 메트릭
        """
        if not target_antigen:
            mapping = self.target_mapping["DEFAULT"]
        else:
            # 타겟명 정규화
            normalized = target_antigen.upper().replace(" ", "").replace("-", "")

            # 매핑 검색
            mapping = None
            for key, value in self.target_mapping.items():
                if key.upper().replace("-", "") in normalized or normalized in key.upper().replace("-", ""):
                    mapping = value
                    break

            if not mapping:
                mapping = self.target_mapping["DEFAULT"]

        gold_standard_key = mapping["gold_standard"]
        gold_standard_data = self.adc_database.get(gold_standard_key)

        if not gold_standard_data:
            # Fallback to Enhertu
            gold_standard_key = "enhertu"
            gold_standard_data = self.adc_database["enhertu"]

        return {
            "key": gold_standard_key,
            "data": gold_standard_data,
            "key_metrics": mapping["key_metrics"],
            "description": mapping["description"]
        }

    def _calculate_relative_score(
        self,
        metric_key: str,
        our_value: float,
        benchmark_value: float
    ) -> RelativeScore:
        """
        상대적 성능 지수 (%) 계산

        Relative Score = (Candidate Value / Benchmark Value) × 100
        (지표 성격에 따라 역산 적용)
        """
        if our_value is None or benchmark_value is None or benchmark_value == 0:
            return RelativeScore(
                metric_name=metric_key,
                our_value=our_value or 0,
                benchmark_value=benchmark_value or 0,
                relative_percent=0,
                direction="unknown",
                interpretation="데이터 부족으로 비교 불가",
                is_favorable=False
            )

        direction = self.metric_direction.get(metric_key, "higher_better")

        # 상대 점수 계산
        if direction == "higher_better":
            relative_percent = (our_value / benchmark_value) * 100
            is_favorable = relative_percent >= 100
        elif direction == "lower_better":
            # 낮을수록 좋은 경우: 역수 적용
            relative_percent = (benchmark_value / our_value) * 100 if our_value > 0 else 0
            is_favorable = our_value <= benchmark_value
        elif direction == "optimal_range":
            # LogP 등 최적 범위가 있는 경우 (2-4 최적)
            optimal_center = 3.0
            our_distance = abs(our_value - optimal_center)
            bench_distance = abs(benchmark_value - optimal_center)
            relative_percent = ((1 - our_distance/optimal_center) / (1 - bench_distance/optimal_center)) * 100 if bench_distance < optimal_center else 100
            is_favorable = our_distance <= bench_distance
        else:
            relative_percent = (our_value / benchmark_value) * 100 if benchmark_value else 0
            is_favorable = True

        # 해석 문장 생성
        interpretation = self._generate_metric_interpretation(
            metric_key, our_value, benchmark_value, relative_percent, is_favorable
        )

        return RelativeScore(
            metric_name=metric_key,
            our_value=our_value,
            benchmark_value=benchmark_value,
            relative_percent=round(relative_percent, 1),
            direction=direction,
            interpretation=interpretation,
            is_favorable=is_favorable
        )

    def _generate_metric_interpretation(
        self,
        metric_key: str,
        our_value: float,
        benchmark_value: float,
        relative_percent: float,
        is_favorable: bool
    ) -> str:
        """지표별 상세 해석 문장 생성"""
        metric_names = {
            "plasma_stability_days": "혈장 안정성",
            "therapeutic_index": "치료 지수",
            "sa_score": "합성 접근성",
            "logp_payload": "LogP",
            "mw_kda": "분자량",
            "dar": "DAR",
            "herg_risk": "hERG 리스크"
        }
        metric_name = metric_names.get(metric_key, metric_key)

        if metric_key == "plasma_stability_days":
            diff = our_value - benchmark_value
            if diff > 0:
                return f"{metric_name}이 기준 약물 대비 {diff:.1f}일 증가하여 순환 반감기가 개선됨. 투여 빈도 감소 가능성 시사."
            elif diff < 0:
                return f"{metric_name}이 기준 약물 대비 {abs(diff):.1f}일 감소. 링커 안정성 보완 필요."
            else:
                return f"{metric_name}이 기준 약물과 동등 수준."

        elif metric_key == "therapeutic_index":
            if relative_percent >= 100:
                return f"치료 지수가 기준 약물의 {relative_percent:.0f}% 수준으로, 더 넓은 안전역(safety margin) 확보."
            else:
                return f"치료 지수가 기준 약물의 {relative_percent:.0f}% 수준. 용량 최적화 필요."

        elif metric_key == "sa_score":
            if our_value < benchmark_value:
                return f"합성 접근성 점수 {our_value:.1f}로, 기준 약물({benchmark_value:.1f}) 대비 합성 용이. 대량 생산 시 비용 절감 기대."
            else:
                return f"합성 접근성 점수 {our_value:.1f}로, 기준 약물({benchmark_value:.1f}) 대비 다소 복잡. 공정 최적화 고려."

        elif metric_key == "logp_payload":
            if 2 <= our_value <= 4:
                return f"페이로드 LogP {our_value:.1f}는 최적 범위(2-4) 내로, 막 투과성과 용해도 균형 양호."
            elif our_value > 4:
                return f"페이로드 LogP {our_value:.1f}로 다소 높음. 조직 축적 가능성 모니터링 필요."
            else:
                return f"페이로드 LogP {our_value:.1f}로 친수성 경향. 세포 침투 효율 검토 필요."

        elif metric_key == "mw_kda":
            diff_percent = ((our_value - benchmark_value) / benchmark_value) * 100
            if abs(diff_percent) < 5:
                return f"분자량이 기준 약물과 유사한 수준({our_value:.0f} vs {benchmark_value:.0f} kDa)."
            elif diff_percent > 0:
                return f"분자량이 기준 약물 대비 {diff_percent:.0f}% 증가하여 조직 침투력은 {100-diff_percent/2:.0f}% 수준으로 예측."
            else:
                return f"분자량이 기준 약물 대비 {abs(diff_percent):.0f}% 감소. 조직 침투력 개선 기대."

        else:
            if is_favorable:
                return f"{metric_name}이 기준 약물 대비 우수함 ({relative_percent:.0f}%)."
            else:
                return f"{metric_name}이 기준 약물 대비 {relative_percent:.0f}% 수준."

    async def analyze(
        self,
        our_adc_metrics: Dict[str, Any],
        target_antigen: str = None,
        payload_class: str = None
    ) -> Dict[str, Any]:
        """
        우리 ADC와 기존 승인 ADC 비교 분석 (동적 벤치마킹)

        Args:
            our_adc_metrics: 우리 ADC의 메트릭
            target_antigen: 타겟 항원
            payload_class: 페이로드 클래스

        Returns:
            벤치마크 비교 결과 + 상대적 성능 지수
        """
        # 1. Gold Standard 자동 선정
        gold_standard = self._get_gold_standard(target_antigen)
        gold_standard_metrics = gold_standard["data"]["metrics"]

        # 2. 추가 비교 대상 선정
        comparators = self._select_comparators(target_antigen, payload_class)
        benchmark_results = []
        our_advantages = []
        our_disadvantages = []

        comparison_properties = [
            ("dar", "DAR", "context_dependent", False),
            ("plasma_stability_days", "Plasma Stability", "higher_better", True),
            ("logp_payload", "LogP (Payload)", "lower_better", True),
            ("sa_score", "Synthetic Accessibility", "lower_better", True),
            ("herg_risk", "hERG Risk", "lower_better", True),
            ("therapeutic_index", "Therapeutic Index", "higher_better", True),
        ]

        for prop_key, prop_name, direction, matters in comparison_properties:
            result = self._compare_property(
                prop_key, prop_name, direction,
                our_adc_metrics, comparators
            )
            benchmark_results.append(result)

            if matters and result.winner == "our_adc":
                our_advantages.append(f"{prop_name}: {result.winner_margin}")
            elif matters and result.winner not in ("our_adc", "tie", "unknown"):
                our_disadvantages.append(f"{prop_name}: vs {result.winner}")

        # 3. 상대적 성능 지수 계산 (Gold Standard 기준)
        relative_scores = {}
        key_metrics = ["plasma_stability_days", "therapeutic_index", "sa_score", "logp_payload", "mw_kda"]

        for metric_key in key_metrics:
            our_value = our_adc_metrics.get(metric_key)
            benchmark_value = gold_standard_metrics.get(metric_key)
            if our_value is not None and benchmark_value is not None:
                relative_scores[metric_key] = self._calculate_relative_score(
                    metric_key, our_value, benchmark_value
                )

        # 4. 비교 테이블 구성
        comparison_table = self._build_comparison_table(our_adc_metrics, comparators)
        winner_analysis = [
            {
                "property": r.property_name,
                "winner": r.winner,
                "margin": r.winner_margin,
                "interpretation": r.interpretation
            }
            for r in benchmark_results
        ]

        # 5. 종합 분석
        competitive_summary = self._generate_competitive_summary(
            our_advantages, our_disadvantages, comparators
        )

        # 6. 상대적 성능 기반 정밀 논리 생성
        relative_reasoning = self._generate_relative_reasoning(
            our_adc_metrics, gold_standard, relative_scores
        )
        reasoning_logic = self._generate_reasoning_logic(
            our_adc_metrics, comparators, our_advantages, our_disadvantages
        )

        # 전체 reasoning 통합
        full_reasoning = f"{relative_reasoning}\n\n{reasoning_logic}"

        confidence_score = self._calculate_confidence(comparators, our_adc_metrics)

        # 7. 레이더 차트용 데이터 준비
        radar_chart_data = self._prepare_radar_chart_data(our_adc_metrics, gold_standard)

        return {
            "success": True,
            # Gold Standard 정보
            "gold_standard": {
                "name": gold_standard["key"],
                "trade_name": gold_standard["data"]["trade_name"],
                "target": gold_standard["data"]["target"],
                "description": gold_standard["description"],
                "metrics": gold_standard_metrics
            },
            # 상대적 성능 지수 (핵심 추가)
            "relative_scores": {
                k: {
                    "our_value": v.our_value,
                    "benchmark_value": v.benchmark_value,
                    "relative_percent": v.relative_percent,
                    "is_favorable": v.is_favorable,
                    "interpretation": v.interpretation
                }
                for k, v in relative_scores.items()
            },
            # 레이더 차트 데이터
            "radar_chart_data": radar_chart_data,
            # 기존 벤치마크 비교
            "benchmark_comparison": {
                "our_adc": our_adc_metrics,
                "competitors": {name: data["metrics"] for name, data in comparators.items()},
                "comparison_table": comparison_table,
                "winner_analysis": winner_analysis
            },
            "competitive_advantages": our_advantages,
            "competitive_disadvantages": our_disadvantages,
            "competitive_summary": competitive_summary,
            "comparator_details": {
                name: {
                    "trade_name": data["trade_name"],
                    "manufacturer": data["manufacturer"],
                    "target": data["target"],
                    "approval_year": data["approval_year"]
                }
                for name, data in comparators.items()
            },
            "reasoning_logic": full_reasoning,
            "confidence_score": confidence_score,
            "confidence_reasoning": f"Gold Standard ({gold_standard['data']['trade_name']}) 기준 + {len(comparators)}개 승인 ADC 비교",
            "referenced_pmids": self._get_reference_pmids(list(comparators.keys()) + [gold_standard["key"]])
        }

    def _generate_relative_reasoning(
        self,
        our_metrics: Dict,
        gold_standard: Dict,
        relative_scores: Dict[str, RelativeScore]
    ) -> str:
        """상대적 성능 기반 정밀 논리 생성"""
        gs_name = gold_standard["data"]["trade_name"]
        parts = [f"[Gold Standard: {gs_name} 대비 상대적 성능 분석]"]

        favorable_count = sum(1 for rs in relative_scores.values() if rs.is_favorable)
        total_count = len(relative_scores)

        parts.append(f"\n본 후보물질은 {gs_name} 대비 {favorable_count}/{total_count}개 핵심 지표에서 우수하거나 동등한 성능을 보입니다.\n")

        for metric_key, rs in relative_scores.items():
            parts.append(f"• {rs.interpretation}")

        # 종합 결론
        if favorable_count >= total_count * 0.7:
            parts.append(f"\n[결론] 본 ADC 후보물질은 {gs_name} 대비 전반적으로 우수한 프로파일을 보이며, 임상 개발 진행에 긍정적입니다.")
        elif favorable_count >= total_count * 0.5:
            parts.append(f"\n[결론] 본 ADC 후보물질은 {gs_name}과 유사한 수준의 효능이 예상되며, 일부 지표에서 차별화 가능성이 있습니다.")
        else:
            parts.append(f"\n[결론] 일부 지표에서 {gs_name} 대비 개선이 필요합니다. 구조 최적화를 통한 보완이 권고됩니다.")

        return "\n".join(parts)

    def _prepare_radar_chart_data(
        self,
        our_metrics: Dict,
        gold_standard: Dict
    ) -> Dict:
        """레이더 차트용 데이터 준비 (100% = Gold Standard)"""
        gs_metrics = gold_standard["data"]["metrics"]

        # 비교할 지표들 (정규화하여 0-150% 스케일로)
        chart_metrics = [
            ("plasma_stability_days", "혈장 안정성", "higher_better"),
            ("therapeutic_index", "치료 지수", "higher_better"),
            ("sa_score", "합성 용이성", "lower_better"),  # 역수
            ("logp_payload", "막 투과성", "optimal_range"),
        ]

        our_normalized = []
        gs_normalized = []  # 항상 100%
        labels = []

        for metric_key, label, direction in chart_metrics:
            our_val = our_metrics.get(metric_key, 0)
            gs_val = gs_metrics.get(metric_key, 1)

            if gs_val == 0:
                gs_val = 1

            if direction == "higher_better":
                our_pct = (our_val / gs_val) * 100
            elif direction == "lower_better":
                our_pct = (gs_val / our_val) * 100 if our_val > 0 else 100
            else:
                # optimal_range: LogP 최적 3.0
                optimal = 3.0
                our_dist = abs(our_val - optimal)
                gs_dist = abs(gs_val - optimal)
                our_pct = ((optimal - our_dist) / (optimal - gs_dist)) * 100 if gs_dist < optimal else 100

            our_normalized.append(min(max(our_pct, 0), 150))  # 0-150% 범위 제한
            gs_normalized.append(100)
            labels.append(label)

        return {
            "labels": labels,
            "our_adc": our_normalized,
            "gold_standard": gs_normalized,
            "gold_standard_name": gold_standard["data"]["trade_name"]
        }

    def _select_comparators(
        self,
        target_antigen: str = None,
        payload_class: str = None,
        max_count: int = 3
    ) -> Dict[str, Dict]:
        """비교 대상 ADC 선정"""
        selected = {}

        if target_antigen:
            for name, data in self.adc_database.items():
                if target_antigen.upper() in data["target"].upper():
                    selected[name] = data
                    if len(selected) >= max_count:
                        return selected

        if payload_class and len(selected) < max_count:
            for name, data in self.adc_database.items():
                if name not in selected:
                    if payload_class.lower() in data["metrics"]["payload_class"].lower():
                        selected[name] = data
                        if len(selected) >= max_count:
                            return selected

        priority_adcs = ["enhertu", "kadcyla", "adcetris"]
        for name in priority_adcs:
            if name not in selected and name in self.adc_database:
                selected[name] = self.adc_database[name]
                if len(selected) >= max_count:
                    return selected

        return selected

    def _compare_property(
        self,
        prop_key: str,
        prop_name: str,
        direction: str,
        our_metrics: Dict,
        comparators: Dict
    ) -> BenchmarkResult:
        """단일 속성 비교"""
        our_value = our_metrics.get(prop_key)
        comp_values = {
            name: data["metrics"].get(prop_key)
            for name, data in comparators.items()
        }

        if prop_key == "herg_risk":
            risk_order = {"Low": 1, "Medium": 2, "High": 3}
            our_numeric = risk_order.get(our_value, 2)
            comp_numeric = {
                name: risk_order.get(val, 2)
                for name, val in comp_values.items()
            }

            best_comp_name = min(comp_numeric, key=comp_numeric.get)
            best_comp_value = comp_numeric[best_comp_name]

            if our_numeric < best_comp_value:
                winner = "our_adc"
                margin = "Lower risk category"
            elif our_numeric == best_comp_value:
                winner = "tie"
                margin = "Same risk level"
            else:
                winner = best_comp_name
                margin = "Higher risk category"

            interpretation = f"Our ADC: {our_value}, Best competitor: {comp_values[best_comp_name]}"

        else:
            if our_value is None:
                return BenchmarkResult(
                    property_name=prop_name,
                    our_value="N/A",
                    competitor_values=comp_values,
                    winner="unknown",
                    winner_margin="Data unavailable",
                    interpretation="Insufficient data for comparison"
                )

            if direction == "higher_better":
                best_comp_name = max(comp_values, key=lambda k: comp_values[k] or 0)
                best_comp_value = comp_values[best_comp_name]
                is_better = our_value > best_comp_value if best_comp_value else True
                margin_value = our_value - best_comp_value if best_comp_value else our_value
            elif direction == "lower_better":
                best_comp_name = min(comp_values, key=lambda k: comp_values[k] or float('inf'))
                best_comp_value = comp_values[best_comp_name]
                is_better = our_value < best_comp_value if best_comp_value else True
                margin_value = best_comp_value - our_value if best_comp_value else 0
            else:
                best_comp_name = list(comp_values.keys())[0]
                best_comp_value = comp_values[best_comp_name]
                is_better = True
                margin_value = 0

            if is_better:
                winner = "our_adc"
                margin = f"+{abs(margin_value):.1f}" if margin_value else "Better"
            elif our_value == best_comp_value:
                winner = "tie"
                margin = "Equivalent"
            else:
                winner = best_comp_name
                margin = f"-{abs(margin_value):.1f}"

            interpretation = f"Our ADC: {our_value}, Best competitor ({best_comp_name}): {best_comp_value}"

        return BenchmarkResult(
            property_name=prop_name,
            our_value=our_value,
            competitor_values=comp_values,
            winner=winner,
            winner_margin=margin,
            interpretation=interpretation
        )

    def _build_comparison_table(
        self,
        our_metrics: Dict,
        comparators: Dict
    ) -> List[Dict]:
        """비교 테이블 구성"""
        properties = [
            ("dar", "DAR"),
            ("mw_kda", "MW (kDa)"),
            ("plasma_stability_days", "Plasma Stability (days)"),
            ("logp_payload", "LogP (Payload)"),
            ("sa_score", "SA Score"),
            ("herg_risk", "hERG Risk"),
            ("therapeutic_index", "Therapeutic Index"),
            ("linker_type", "Linker Type"),
            ("payload_class", "Payload Class")
        ]

        table = []
        for prop_key, prop_name in properties:
            row = {
                "property": prop_name,
                "our_adc": our_metrics.get(prop_key, "N/A")
            }
            for comp_name, comp_data in comparators.items():
                row[comp_name] = comp_data["metrics"].get(prop_key, "N/A")
            table.append(row)

        return table

    def _generate_competitive_summary(
        self,
        advantages: List[str],
        disadvantages: List[str],
        comparators: Dict
    ) -> str:
        """경쟁 요약 생성"""
        comp_names = [data["trade_name"] for data in comparators.values()]
        summary_parts = []

        if advantages:
            summary_parts.append(
                f"Our ADC demonstrates superiority in {len(advantages)} key metrics "
                f"compared to {', '.join(comp_names)}: {'; '.join(advantages[:3])}."
            )

        if disadvantages:
            summary_parts.append(
                f"Areas requiring attention: {'; '.join(disadvantages[:2])}."
            )

        if not advantages and not disadvantages:
            summary_parts.append(
                f"Our ADC shows comparable performance to approved ADCs ({', '.join(comp_names)})."
            )

        return " ".join(summary_parts)

    def _generate_reasoning_logic(
        self,
        our_metrics: Dict,
        comparators: Dict,
        advantages: List[str],
        disadvantages: List[str]
    ) -> str:
        """추론 로직 생성"""
        parts = []

        if our_metrics.get("plasma_stability_days", 0) > 6:
            parts.append(
                f"혈장 안정성 {our_metrics.get('plasma_stability_days')}일은 "
                "ADC의 순환 반감기를 개선하여 투여 빈도를 줄일 수 있는 이점을 제공함."
            )

        sa_score = our_metrics.get("sa_score", 5)
        if sa_score < 4:
            parts.append(
                f"합성 접근성 점수 {sa_score}는 표준 유기합성 실험실에서 "
                "합성 가능한 수준으로, 대량 생산 시 비용 절감이 기대됨."
            )

        herg = our_metrics.get("herg_risk", "Medium")
        if herg == "Low":
            parts.append(
                "hERG 리스크가 'Low'로 분류되어 심장 독성 우려가 낮음."
            )

        comp_count = len(comparators)
        adv_count = len(advantages)
        parts.append(
            f"기존 승인 ADC {comp_count}개와의 비교 결과, "
            f"{adv_count}개 핵심 지표에서 우위를 보임."
        )

        return " ".join(parts)

    def _calculate_confidence(self, comparators: Dict, our_metrics: Dict) -> float:
        """신뢰도 점수 계산"""
        base_score = 70
        base_score += len(comparators) * 5

        required_metrics = ["dar", "plasma_stability_days", "sa_score", "herg_risk"]
        present_count = sum(1 for m in required_metrics if our_metrics.get(m) is not None)
        base_score += (present_count / len(required_metrics)) * 10

        return min(100, base_score)

    def _get_reference_pmids(self, comparator_names: List[str]) -> List[str]:
        """참조 PMID 목록"""
        pmid_database = {
            "enhertu": ["31825569", "33556277", "35172035"],
            "kadcyla": ["22149875", "25085099"],
            "adcetris": ["21296805", "25559804"],
            "padcev": ["31586450", "34530666"],
            "trodelvy": ["31116389", "33882206"]
        }

        pmids = []
        for name in comparator_names:
            pmids.extend(pmid_database.get(name, []))

        return list(set(pmids))[:10]


async def run_benchmark_comparison(
    our_adc_metrics: Dict[str, Any],
    target_antigen: str = None,
    payload_class: str = None
) -> Dict[str, Any]:
    """벤치마크 비교 실행"""
    agent = BenchmarkComparisonAgent()
    return await agent.analyze(our_adc_metrics, target_antigen, payload_class)
