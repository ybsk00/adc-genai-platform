"""
Confidence Calculator - 에이전트 추론 신뢰도 계산
각 섹션의 Confidence Score (0-100%) 산정
"""
from typing import Dict, Any, List, Optional
from dataclasses import dataclass
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class ConfidenceLevel(Enum):
    """신뢰도 레벨"""
    HIGH = "High"      # 80-100%
    MEDIUM = "Medium"  # 60-79%
    LOW = "Low"        # 0-59%


@dataclass
class ConfidenceResult:
    """신뢰도 계산 결과"""
    score: float
    level: ConfidenceLevel
    reasoning: str
    limitations: List[str]
    breakdown: Dict[str, float]


class ConfidenceWeights:
    """가중치 설정"""
    LITERATURE_EVIDENCE = 0.30   # 논문 근거 (30%)
    GOLDEN_SET_MATCH = 0.25      # Golden Set 매칭 (25%)
    PREDICTION_ACCURACY = 0.25  # 예측 모델 정확도 (25%)
    DATA_COMPLETENESS = 0.20    # 데이터 완전성 (20%)


def calculate_confidence(
    pmid_count: int = 0,
    golden_set_matches: int = 0,
    prediction_model_accuracy: float = 0.8,
    data_completeness: float = 1.0,
    custom_factors: Dict[str, float] = None
) -> ConfidenceResult:
    """
    Confidence Score 계산

    Args:
        pmid_count: 참조된 PMID(논문) 개수
        golden_set_matches: Golden Set에서 유사 구조 매칭 개수
        prediction_model_accuracy: 예측 모델 정확도 (0-1)
        data_completeness: 입력 데이터 완전성 (0-1)
        custom_factors: 커스텀 가중 요소

    Returns:
        ConfidenceResult 객체
    """
    weights = ConfidenceWeights()

    # 각 요소별 점수 계산 (0-100 스케일)

    # 1. 논문 근거 점수 (최대 5개 논문 = 100점)
    literature_score = min(pmid_count / 5, 1.0) * 100

    # 2. Golden Set 매칭 점수 (최대 3개 매칭 = 100점)
    golden_set_score = min(golden_set_matches / 3, 1.0) * 100

    # 3. 예측 모델 정확도 점수
    prediction_score = prediction_model_accuracy * 100

    # 4. 데이터 완전성 점수
    completeness_score = data_completeness * 100

    # 가중 평균 계산
    breakdown = {
        "literature_evidence": literature_score,
        "golden_set_match": golden_set_score,
        "prediction_accuracy": prediction_score,
        "data_completeness": completeness_score
    }

    total_score = (
        literature_score * weights.LITERATURE_EVIDENCE +
        golden_set_score * weights.GOLDEN_SET_MATCH +
        prediction_score * weights.PREDICTION_ACCURACY +
        completeness_score * weights.DATA_COMPLETENESS
    )

    # 커스텀 요소 적용
    if custom_factors:
        for factor_name, factor_value in custom_factors.items():
            # 커스텀 요소는 ±10% 범위 내에서 조정
            adjustment = (factor_value - 0.5) * 20  # -10 ~ +10
            total_score += adjustment
            breakdown[factor_name] = factor_value * 100

    # 점수 범위 제한 (0-100)
    total_score = max(0, min(100, total_score))

    # 레벨 결정
    if total_score >= 80:
        level = ConfidenceLevel.HIGH
    elif total_score >= 60:
        level = ConfidenceLevel.MEDIUM
    else:
        level = ConfidenceLevel.LOW

    # 한계점 목록 생성
    limitations = []
    if pmid_count < 2:
        limitations.append("근거 논문 부족 (2개 미만)")
    if golden_set_matches == 0:
        limitations.append("Golden Set에 유사 구조 없음")
    if prediction_model_accuracy < 0.7:
        limitations.append("예측 모델 정확도 제한적 (70% 미만)")
    if data_completeness < 0.8:
        limitations.append("입력 데이터 불완전 (일부 필드 누락)")

    # 추론 근거 생성
    reasoning_parts = []
    if pmid_count > 0:
        reasoning_parts.append(f"{pmid_count}개의 peer-reviewed 논문 참조")
    if golden_set_matches > 0:
        reasoning_parts.append(f"Golden Set에서 {golden_set_matches}개 유사 구조 발견")
    if prediction_model_accuracy >= 0.8:
        reasoning_parts.append(f"예측 모델 정확도 {prediction_model_accuracy*100:.0f}%")
    if data_completeness >= 0.9:
        reasoning_parts.append("완전한 입력 데이터")

    reasoning = "; ".join(reasoning_parts) if reasoning_parts else "제한된 근거 데이터로 추론됨"

    return ConfidenceResult(
        score=round(total_score, 1),
        level=level,
        reasoning=reasoning,
        limitations=limitations,
        breakdown={k: round(v, 1) for k, v in breakdown.items()}
    )


def calculate_section_confidence(
    section_type: str,
    agent_data: Dict[str, Any]
) -> ConfidenceResult:
    """
    보고서 섹션별 Confidence Score 계산

    Args:
        section_type: 섹션 타입 (target_analysis, molecular_design, etc.)
        agent_data: 에이전트 출력 데이터

    Returns:
        ConfidenceResult 객체
    """
    # 에이전트 데이터에서 기본 값 추출
    pmids = agent_data.get('referenced_pmids', [])
    golden_matches = agent_data.get('golden_set_matches', [])

    # 섹션별 특화된 계산
    custom_factors = {}

    if section_type == "target_analysis":
        # 타겟 분석: UniProt 데이터 존재 여부
        if agent_data.get('uniprot_data'):
            custom_factors['uniprot_verified'] = 0.9
        if agent_data.get('binding_affinity'):
            custom_factors['binding_data'] = 0.8

    elif section_type == "molecular_design":
        # 분자 설계: SMILES 유효성, 구조적 유사도
        if agent_data.get('smiles_validated'):
            custom_factors['structure_valid'] = 0.9
        if agent_data.get('structural_similarity', 0) > 0.7:
            custom_factors['similar_approved'] = 0.85

    elif section_type == "synthesis_path":
        # 합성 경로: SA Score, 합성 단계 수
        sa_score = agent_data.get('sa_score', 5)
        custom_factors['synthetic_feasibility'] = max(0, 1 - sa_score / 10)

    elif section_type == "toxicology":
        # 독성 분석: hERG, 간독성 예측
        toxicity_data = agent_data.get('toxicity_predictions', {})
        if toxicity_data:
            herg_risk = toxicity_data.get('herg_risk', 'Medium')
            custom_factors['safety_profile'] = {
                'Low': 0.9, 'Medium': 0.6, 'High': 0.3
            }.get(herg_risk, 0.5)

    elif section_type == "benchmark_comparison":
        # 벤치마크 비교: 비교 데이터 완전성
        competitors = agent_data.get('competitors', {})
        if len(competitors) >= 3:
            custom_factors['comprehensive_comparison'] = 0.9
        elif len(competitors) >= 1:
            custom_factors['limited_comparison'] = 0.7

    elif section_type == "patent_landscape":
        # 특허 분석: 특허 검색 커버리지
        patents_found = agent_data.get('patents_analyzed', 0)
        custom_factors['patent_coverage'] = min(patents_found / 10, 1.0)

    # 데이터 완전성 계산
    required_fields = _get_required_fields(section_type)
    present_fields = sum(1 for f in required_fields if agent_data.get(f) is not None)
    data_completeness = present_fields / len(required_fields) if required_fields else 1.0

    # 모델 정확도 (에이전트에서 제공하거나 기본값 사용)
    model_accuracy = agent_data.get('model_accuracy', 0.8)

    return calculate_confidence(
        pmid_count=len(pmids),
        golden_set_matches=len(golden_matches),
        prediction_model_accuracy=model_accuracy,
        data_completeness=data_completeness,
        custom_factors=custom_factors
    )


def _get_required_fields(section_type: str) -> List[str]:
    """섹션별 필수 필드 목록"""
    required_fields = {
        "target_analysis": ["target_name", "uniprot_id", "biological_function", "expression_profile"],
        "molecular_design": ["smiles", "linker_type", "payload_class", "conjugation_site"],
        "synthesis_path": ["sa_score", "synthesis_steps", "estimated_yield"],
        "physicochemical": ["molecular_weight", "logp", "hbd", "hba", "tpsa"],
        "toxicology": ["herg_risk", "hepatotoxicity_risk", "off_target_binding"],
        "benchmark_comparison": ["our_metrics", "competitor_metrics", "advantages"],
        "patent_landscape": ["relevant_patents", "freedom_to_operate", "novelty_assessment"],
    }
    return required_fields.get(section_type, [])


def aggregate_report_confidence(section_confidences: Dict[str, ConfidenceResult]) -> ConfidenceResult:
    """
    전체 보고서 종합 Confidence Score

    Args:
        section_confidences: 각 섹션별 ConfidenceResult

    Returns:
        종합 ConfidenceResult
    """
    if not section_confidences:
        return ConfidenceResult(
            score=0,
            level=ConfidenceLevel.LOW,
            reasoning="No section data available",
            limitations=["No analysis data"],
            breakdown={}
        )

    # 섹션별 가중치 (중요도 순)
    section_weights = {
        "molecular_design": 0.25,
        "target_analysis": 0.20,
        "toxicology": 0.20,
        "benchmark_comparison": 0.15,
        "synthesis_path": 0.10,
        "patent_landscape": 0.10
    }

    total_weight = 0
    weighted_score = 0
    all_limitations = []
    breakdown = {}

    for section_name, conf_result in section_confidences.items():
        weight = section_weights.get(section_name, 0.1)
        weighted_score += conf_result.score * weight
        total_weight += weight
        all_limitations.extend(conf_result.limitations)
        breakdown[section_name] = conf_result.score

    # 정규화
    if total_weight > 0:
        final_score = weighted_score / total_weight
    else:
        final_score = sum(c.score for c in section_confidences.values()) / len(section_confidences)

    # 중복 제거된 한계점
    unique_limitations = list(set(all_limitations))

    # 레벨 결정
    if final_score >= 80:
        level = ConfidenceLevel.HIGH
    elif final_score >= 60:
        level = ConfidenceLevel.MEDIUM
    else:
        level = ConfidenceLevel.LOW

    # 종합 추론 근거
    high_conf_sections = [name for name, conf in section_confidences.items() if conf.level == ConfidenceLevel.HIGH]
    low_conf_sections = [name for name, conf in section_confidences.items() if conf.level == ConfidenceLevel.LOW]

    reasoning_parts = []
    if high_conf_sections:
        reasoning_parts.append(f"High confidence in: {', '.join(high_conf_sections)}")
    if low_conf_sections:
        reasoning_parts.append(f"Limited confidence in: {', '.join(low_conf_sections)}")

    reasoning = "; ".join(reasoning_parts) if reasoning_parts else "Mixed confidence across sections"

    return ConfidenceResult(
        score=round(final_score, 1),
        level=level,
        reasoning=reasoning,
        limitations=unique_limitations[:5],  # 상위 5개만
        breakdown=breakdown
    )


class ConfidenceCalculator:
    """Confidence calculation service class"""

    @staticmethod
    def calculate(
        pmid_count: int = 0,
        golden_set_matches: int = 0,
        prediction_model_accuracy: float = 0.8,
        data_completeness: float = 1.0,
        custom_factors: Dict[str, float] = None
    ) -> ConfidenceResult:
        return calculate_confidence(
            pmid_count, golden_set_matches,
            prediction_model_accuracy, data_completeness,
            custom_factors
        )

    @staticmethod
    def calculate_section(section_type: str, agent_data: Dict[str, Any]) -> ConfidenceResult:
        return calculate_section_confidence(section_type, agent_data)

    @staticmethod
    def aggregate(section_confidences: Dict[str, ConfidenceResult]) -> ConfidenceResult:
        return aggregate_report_confidence(section_confidences)
