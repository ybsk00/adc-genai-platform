"""
Visualization Service - 보고서용 자동 차트 생성
Lipinski Chart, DAR Distribution, Linker Cleavage Diagram
"""
import io
import base64
from typing import Dict, Any, List, Optional
import logging

logger = logging.getLogger(__name__)

# Lazy import for matplotlib to avoid startup overhead
_plt = None
_np = None


def _get_matplotlib():
    global _plt
    if _plt is None:
        import matplotlib
        matplotlib.use('Agg')  # Non-GUI backend
        import matplotlib.pyplot as plt
        _plt = plt
    return _plt


def _get_numpy():
    global _np
    if _np is None:
        import numpy as np
        _np = np
    return _np


def fig_to_base64(fig) -> str:
    """Matplotlib figure를 base64 PNG로 변환"""
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    buffer.seek(0)
    img_base64 = base64.b64encode(buffer.getvalue()).decode('utf-8')
    _get_matplotlib().close(fig)
    return img_base64


def generate_lipinski_chart(metrics: Dict[str, float]) -> str:
    """
    Lipinski's Rule of Five 레이더 차트 생성

    Args:
        metrics: {mw, logp, hbd, hba, rotatable_bonds}

    Returns:
        Base64 encoded PNG image
    """
    plt = _get_matplotlib()
    np = _get_numpy()

    categories = ['MW\n(<500)', 'LogP\n(<5)', 'HBD\n(<5)', 'HBA\n(<10)', 'RotBonds\n(<10)']

    # Normalize values (1.0 = at limit, <1.0 = within rule, >1.0 = exceeds rule)
    values = [
        metrics.get('mw', 400) / 500,
        metrics.get('logp', 3) / 5,
        metrics.get('hbd', 2) / 5,
        metrics.get('hba', 5) / 10,
        metrics.get('rotatable_bonds', 5) / 10
    ]

    # Clamp values for display
    display_values = [min(v, 1.5) for v in values]

    # Create radar chart
    angles = np.linspace(0, 2 * np.pi, len(categories), endpoint=False).tolist()
    display_values_closed = display_values + display_values[:1]
    angles_closed = angles + angles[:1]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

    # Draw threshold circle at 1.0
    threshold_circle = [1.0] * (len(categories) + 1)
    ax.plot(angles_closed, threshold_circle, 'r--', linewidth=2, label='Lipinski Limit', alpha=0.7)
    ax.fill(angles_closed, threshold_circle, 'red', alpha=0.05)

    # Determine colors based on compliance
    colors = []
    for v in values:
        if v <= 1.0:
            colors.append('#22c55e')  # Green - compliant
        else:
            colors.append('#ef4444')  # Red - violation

    # Plot values
    ax.fill(angles_closed, display_values_closed, alpha=0.25, color='#3b82f6')
    ax.plot(angles_closed, display_values_closed, 'o-', linewidth=2, color='#3b82f6', markersize=8)

    # Add colored markers for compliance status
    for angle, val, color in zip(angles, display_values, colors):
        ax.plot(angle, val, 'o', color=color, markersize=12, zorder=5)

    # Customize
    ax.set_xticks(angles)
    ax.set_xticklabels(categories, fontsize=10, fontweight='bold')
    ax.set_ylim(0, 1.5)
    ax.set_yticks([0.5, 1.0, 1.5])
    ax.set_yticklabels(['50%', '100%', '150%'], fontsize=8)
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))

    # Title with compliance count
    compliant_count = sum(1 for v in values if v <= 1.0)
    ax.set_title(f"Lipinski's Rule of Five\n({compliant_count}/5 Rules Passed)",
                 fontsize=14, fontweight='bold', pad=20)

    return fig_to_base64(fig)


def generate_dar_distribution(target_dar: int, std_dev: float = 0.8) -> str:
    """
    DAR 분포 예측 그래프 생성

    Args:
        target_dar: 목표 DAR 값
        std_dev: 분포 표준편차

    Returns:
        Base64 encoded PNG image
    """
    plt = _get_matplotlib()
    np = _get_numpy()
    from scipy.stats import norm

    dar_values = np.arange(0, 10, 0.05)
    distribution = norm.pdf(dar_values, target_dar, std_dev)

    fig, ax = plt.subplots(figsize=(8, 4))

    # Fill distribution
    ax.fill_between(dar_values, distribution, alpha=0.3, color='#3b82f6')
    ax.plot(dar_values, distribution, color='#3b82f6', linewidth=2)

    # Mark target DAR
    ax.axvline(target_dar, color='#ef4444', linestyle='--', linewidth=2,
               label=f'Target DAR: {target_dar}')

    # Add shaded regions for different DAR ranges
    ax.axvspan(0, 2, alpha=0.1, color='yellow', label='Low DAR (0-2)')
    ax.axvspan(target_dar - 1, target_dar + 1, alpha=0.2, color='green', label='Optimal Range')
    ax.axvspan(8, 10, alpha=0.1, color='red', label='High DAR (8+)')

    # Calculate probability within optimal range
    prob_optimal = norm.cdf(target_dar + 1, target_dar, std_dev) - norm.cdf(target_dar - 1, target_dar, std_dev)

    ax.set_xlabel('Drug-to-Antibody Ratio (DAR)', fontsize=12)
    ax.set_ylabel('Probability Density', fontsize=12)
    ax.set_title(f'Predicted DAR Distribution\n(Optimal Range Probability: {prob_optimal*100:.1f}%)',
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlim(0, 10)
    ax.grid(True, alpha=0.3)

    return fig_to_base64(fig)


def generate_linker_cleavage_diagram(linker_type: str, smiles: str) -> str:
    """
    링커 절단 부위 시뮬레이션 다이어그램
    RDKit을 사용하여 분자 구조와 절단 부위 표시

    Args:
        linker_type: 링커 타입 (cleavable, non-cleavable)
        smiles: 링커 SMILES

    Returns:
        Base64 encoded PNG image
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
        import re

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return _generate_placeholder_cleavage_diagram(linker_type)

        AllChem.Compute2DCoords(mol)

        # Find potential cleavage sites based on linker type
        highlight_atoms = []
        highlight_bonds = []

        if linker_type.lower() == 'cleavable':
            # Look for peptide bonds (amide), disulfide bridges, hydrazone
            patterns = [
                ('C(=O)N', 'Peptide Bond'),
                ('SS', 'Disulfide'),
                ('C=NN', 'Hydrazone'),
                ('OC(=O)', 'Ester'),
            ]

            for pattern, name in patterns:
                substructure = Chem.MolFromSmarts(pattern)
                if substructure:
                    matches = mol.GetSubstructMatches(substructure)
                    for match in matches:
                        highlight_atoms.extend(match)

        # Create drawing
        drawer = rdMolDraw2D.MolDraw2DCairo(500, 350)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.drawOptions().addAtomIndices = False

        # Highlight cleavage sites in red
        if highlight_atoms:
            highlight_colors = {atom: (1.0, 0.2, 0.2) for atom in highlight_atoms}
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms,
                               highlightAtomColors=highlight_colors)
        else:
            drawer.DrawMolecule(mol)

        drawer.FinishDrawing()
        img_data = drawer.GetDrawingText()

        return base64.b64encode(img_data).decode('utf-8')

    except Exception as e:
        logger.warning(f"RDKit drawing failed: {e}")
        return _generate_placeholder_cleavage_diagram(linker_type)


def _generate_placeholder_cleavage_diagram(linker_type: str) -> str:
    """RDKit 실패 시 대체 다이어그램"""
    plt = _get_matplotlib()

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, f"Linker Type: {linker_type}\n\n[Molecular structure\nrendering unavailable]",
            ha='center', va='center', fontsize=12, transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='#f3f4f6', edgecolor='#d1d5db'))
    ax.axis('off')
    ax.set_title(f'Linker Cleavage Analysis', fontsize=14, fontweight='bold')

    return fig_to_base64(fig)


def generate_benchmark_comparison_chart(our_adc: Dict, competitors: Dict) -> str:
    """
    벤치마크 비교 막대 차트 생성

    Args:
        our_adc: 우리 ADC의 메트릭
        competitors: 경쟁 ADC들의 메트릭

    Returns:
        Base64 encoded PNG image
    """
    plt = _get_matplotlib()
    np = _get_numpy()

    # Properties to compare (normalized to 0-100 scale)
    properties = ['Plasma\nStability', 'Therapeutic\nIndex', 'Synthetic\nAccessibility', 'Safety\nProfile']

    # Normalize scores (higher = better)
    our_scores = [
        min(our_adc.get('plasma_stability_days', 5) / 10 * 100, 100),
        min(our_adc.get('therapeutic_index', 10) / 20 * 100, 100),
        max(0, 100 - our_adc.get('sa_score', 3) * 10),  # Lower SA = better
        100 - ({'High': 80, 'Medium': 50, 'Low': 20}.get(our_adc.get('herg_risk', 'Medium'), 50))
    ]

    competitor_scores = {}
    for name, data in competitors.items():
        competitor_scores[name] = [
            min(data.get('plasma_stability_days', 5) / 10 * 100, 100),
            min(data.get('therapeutic_index', 10) / 20 * 100, 100),
            max(0, 100 - data.get('sa_score', 4) * 10),
            100 - ({'High': 80, 'Medium': 50, 'Low': 20}.get(data.get('herg_risk', 'Medium'), 50))
        ]

    x = np.arange(len(properties))
    width = 0.2

    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot our ADC
    bars1 = ax.bar(x - width * 1.5, our_scores, width, label='Our ADC', color='#3b82f6', edgecolor='white')

    # Plot competitors
    colors = ['#f97316', '#22c55e', '#a855f7']
    for i, (name, scores) in enumerate(competitor_scores.items()):
        offset = x - width * 0.5 + width * i
        ax.bar(offset, scores, width, label=name.title(), color=colors[i % len(colors)], edgecolor='white')

    ax.set_ylabel('Score (0-100)', fontsize=12)
    ax.set_title('Benchmark Comparison vs Approved ADCs', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(properties, fontsize=10)
    ax.legend(loc='upper right', fontsize=9)
    ax.set_ylim(0, 110)
    ax.grid(True, axis='y', alpha=0.3)

    # Add value labels
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.0f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=8)

    return fig_to_base64(fig)


def generate_confidence_gauge(score: float, section_name: str) -> str:
    """
    Confidence Score 게이지 차트 생성

    Args:
        score: 0-100 신뢰도 점수
        section_name: 섹션 이름

    Returns:
        Base64 encoded PNG image
    """
    plt = _get_matplotlib()
    np = _get_numpy()

    fig, ax = plt.subplots(figsize=(4, 2.5))

    # Determine color based on score
    if score >= 80:
        color = '#22c55e'  # Green
        level = 'High'
    elif score >= 60:
        color = '#f59e0b'  # Amber
        level = 'Medium'
    else:
        color = '#ef4444'  # Red
        level = 'Low'

    # Create semi-circle gauge
    theta = np.linspace(0, np.pi, 100)

    # Background arc
    ax.fill_between(theta, 0.7, 1.0, alpha=0.2, color='#e5e7eb')

    # Colored arc based on score
    score_angle = np.pi * score / 100
    filled_theta = np.linspace(0, score_angle, int(score))
    ax.fill_between(filled_theta, 0.7, 1.0, alpha=0.8, color=color)

    # Convert to cartesian for plotting
    ax.set_xlim(0, np.pi)
    ax.set_ylim(0, 1.2)

    # Add score text
    ax.text(np.pi / 2, 0.4, f'{score:.0f}%', ha='center', va='center',
            fontsize=24, fontweight='bold', color=color)
    ax.text(np.pi / 2, 0.15, f'Confidence: {level}', ha='center', va='center',
            fontsize=10, color='#6b7280')

    ax.axis('off')
    ax.set_title(section_name, fontsize=11, fontweight='bold', pad=10)

    return fig_to_base64(fig)


def generate_sa_score_chart(sa_score: float, comparison_data: Dict = None) -> str:
    """
    Synthetic Accessibility Score 비교 차트

    Args:
        sa_score: SA Score (1-10, lower = easier)
        comparison_data: 경쟁 약물들과의 비교 데이터

    Returns:
        Base64 encoded PNG image
    """
    plt = _get_matplotlib()
    np = _get_numpy()

    fig, ax = plt.subplots(figsize=(8, 3))

    # Create horizontal bar showing SA score scale
    categories = ['Very Easy\n(1-2)', 'Easy\n(2-4)', 'Moderate\n(4-6)', 'Difficult\n(6-8)', 'Very Difficult\n(8-10)']
    colors = ['#22c55e', '#84cc16', '#eab308', '#f97316', '#ef4444']

    # Scale positions
    positions = [1.5, 3, 5, 7, 9]

    # Draw scale background
    for i, (pos, color) in enumerate(zip([0, 2, 4, 6, 8], colors)):
        ax.barh(0, 2, left=pos, height=0.3, color=color, alpha=0.3)

    # Mark our SA score
    ax.plot(sa_score, 0, 'v', markersize=20, color='#3b82f6', label=f'Our ADC: {sa_score:.1f}')
    ax.annotate(f'{sa_score:.1f}', (sa_score, 0.25), ha='center', fontsize=12, fontweight='bold', color='#3b82f6')

    # Add comparison markers if provided
    if comparison_data:
        markers = ['o', 's', '^', 'd']
        comp_colors = ['#f97316', '#22c55e', '#a855f7']
        for i, (name, data) in enumerate(comparison_data.items()):
            comp_sa = data.get('sa_score', 5)
            ax.plot(comp_sa, 0, markers[i % 4], markersize=12, color=comp_colors[i % 3],
                   label=f'{name.title()}: {comp_sa:.1f}', alpha=0.7)

    ax.set_xlim(0, 10)
    ax.set_ylim(-0.5, 0.8)
    ax.set_xlabel('Synthetic Accessibility Score', fontsize=12)
    ax.set_title('Synthetic Accessibility Comparison', fontsize=14, fontweight='bold')
    ax.set_yticks([])
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=4, fontsize=9)

    # Add scale labels
    for pos, label in zip(positions, categories):
        ax.text(pos, -0.35, label, ha='center', va='top', fontsize=8)

    plt.tight_layout()
    return fig_to_base64(fig)


def generate_relative_radar_chart(radar_data: Dict) -> str:
    """
    상대적 성능 레이더 차트 생성
    중앙 100% = Gold Standard, 바깥쪽 = 우수, 안쪽 = 열등

    Args:
        radar_data: {
            "labels": ["혈장 안정성", "치료 지수", ...],
            "our_adc": [110, 95, 120, ...],  # % 값
            "gold_standard": [100, 100, 100, ...],
            "gold_standard_name": "Enhertu"
        }

    Returns:
        Base64 encoded PNG image
    """
    plt = _get_matplotlib()
    np = _get_numpy()

    labels = radar_data.get("labels", [])
    our_values = radar_data.get("our_adc", [])
    gs_values = radar_data.get("gold_standard", [])
    gs_name = radar_data.get("gold_standard_name", "Gold Standard")

    if not labels or len(labels) < 3:
        # 최소 3개 지표 필요
        return _generate_placeholder_radar()

    num_vars = len(labels)

    # 각도 계산
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1]  # 닫힌 도형을 위해

    # 데이터 닫기
    our_values_closed = our_values + our_values[:1]
    gs_values_closed = gs_values + gs_values[:1]

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))

    # 100% 기준선 (Gold Standard)
    ax.plot(angles, gs_values_closed, 'o-', linewidth=2, color='#f97316',
            label=f'{gs_name} (100%)', markersize=8)
    ax.fill(angles, gs_values_closed, alpha=0.15, color='#f97316')

    # 우리 ADC
    ax.plot(angles, our_values_closed, 'o-', linewidth=3, color='#3b82f6',
            label='Our ADC', markersize=10)
    ax.fill(angles, our_values_closed, alpha=0.25, color='#3b82f6')

    # 레이블 설정
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=11, fontweight='bold')

    # Y축 설정 (0-150% 범위)
    ax.set_ylim(0, 150)
    ax.set_yticks([50, 100, 150])
    ax.set_yticklabels(['50%', '100%', '150%'], fontsize=9, color='#6b7280')

    # 100% 원 강조
    circle_100 = plt.Circle((0, 0), 100, transform=ax.transData._b, fill=False,
                            color='#ef4444', linewidth=2, linestyle='--', alpha=0.5)
    # Note: polar axes에서는 직접 원을 그리기 어려움, 대신 y=100 gridline 사용

    # 범례
    ax.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize=10)

    # 제목
    ax.set_title(f'Relative Performance vs {gs_name}\n(100% = Gold Standard)',
                 fontsize=14, fontweight='bold', pad=20)

    # 각 점에 % 값 표시
    for angle, val in zip(angles[:-1], our_values):
        # 100% 대비 색상 결정
        if val >= 100:
            color = '#22c55e'  # 녹색 (우수)
            offset = 15
        else:
            color = '#ef4444'  # 빨강 (열등)
            offset = -20

        ax.annotate(f'{val:.0f}%',
                    xy=(angle, val),
                    xytext=(5, offset),
                    textcoords='offset points',
                    ha='center',
                    fontsize=10,
                    fontweight='bold',
                    color=color)

    plt.tight_layout()
    return fig_to_base64(fig)


def _generate_placeholder_radar() -> str:
    """레이더 차트 대체 이미지"""
    plt = _get_matplotlib()

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, "Insufficient data\nfor radar chart",
            ha='center', va='center', fontsize=12, transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='#f3f4f6', edgecolor='#d1d5db'))
    ax.axis('off')
    ax.set_title('Relative Performance', fontsize=14, fontweight='bold')

    return fig_to_base64(fig)


def generate_relative_scores_bar_chart(relative_scores: Dict) -> str:
    """
    상대적 성능 지수 막대 차트
    100%를 기준선으로 하고 각 지표의 상대값 표시

    Args:
        relative_scores: {
            "metric_name": {
                "relative_percent": 115.5,
                "is_favorable": True,
                "interpretation": "..."
            }
        }
    """
    plt = _get_matplotlib()
    np = _get_numpy()

    metrics = list(relative_scores.keys())
    values = [relative_scores[m].get("relative_percent", 100) for m in metrics]
    favorable = [relative_scores[m].get("is_favorable", True) for m in metrics]

    # 한글 지표명 매핑
    metric_labels = {
        "plasma_stability_days": "혈장 안정성",
        "therapeutic_index": "치료 지수",
        "sa_score": "합성 용이성",
        "logp_payload": "막 투과성",
        "mw_kda": "분자량"
    }
    labels = [metric_labels.get(m, m) for m in metrics]

    # 색상 결정 (100% 이상이면 녹색, 미만이면 빨강)
    colors = ['#22c55e' if f else '#ef4444' for f in favorable]

    fig, ax = plt.subplots(figsize=(10, 5))

    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, values, color=colors, edgecolor='white', height=0.6)

    # 100% 기준선
    ax.axvline(x=100, color='#1e3a5f', linewidth=2, linestyle='--', label='Gold Standard (100%)')

    # 값 라벨
    for i, (bar, val, fav) in enumerate(zip(bars, values, favorable)):
        text_color = '#166534' if fav else '#991b1b'
        prefix = '+' if val > 100 else ''
        ax.text(bar.get_width() + 2, bar.get_y() + bar.get_height()/2,
                f'{prefix}{val-100:.1f}%' if val != 100 else '±0%',
                va='center', fontsize=11, fontweight='bold', color=text_color)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel('Gold Standard 대비 상대적 성능 (%)', fontsize=12)
    ax.set_title('Relative Efficiency Score vs Gold Standard', fontsize=14, fontweight='bold')
    ax.set_xlim(0, max(values) + 30)
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, axis='x', alpha=0.3)

    plt.tight_layout()
    return fig_to_base64(fig)


class VisualizationService:
    """Report visualization service"""

    @staticmethod
    def generate_lipinski_chart(metrics: Dict[str, float]) -> str:
        return generate_lipinski_chart(metrics)

    @staticmethod
    def generate_dar_distribution(target_dar: int, std_dev: float = 0.8) -> str:
        return generate_dar_distribution(target_dar, std_dev)

    @staticmethod
    def generate_linker_cleavage_diagram(linker_type: str, smiles: str) -> str:
        return generate_linker_cleavage_diagram(linker_type, smiles)

    @staticmethod
    def generate_benchmark_comparison_chart(our_adc: Dict, competitors: Dict) -> str:
        return generate_benchmark_comparison_chart(our_adc, competitors)

    @staticmethod
    def generate_confidence_gauge(score: float, section_name: str) -> str:
        return generate_confidence_gauge(score, section_name)

    @staticmethod
    def generate_sa_score_chart(sa_score: float, comparison_data: Dict = None) -> str:
        return generate_sa_score_chart(sa_score, comparison_data)

    @staticmethod
    def generate_relative_radar_chart(radar_data: Dict) -> str:
        return generate_relative_radar_chart(radar_data)

    @staticmethod
    def generate_relative_scores_bar_chart(relative_scores: Dict) -> str:
        return generate_relative_scores_bar_chart(relative_scores)
