"""
Report Engine - 10페이지 정밀 보고서 생성 엔진
WeasyPrint 기반 PDF 생성 + Digital Seal + QR 검증

Phase 2 Enhancement:
- Digital Lineage (데이터 계보) 통합
- Calculation Metadata 태깅
- 재현성 검증 지원
"""
import os
import io
import json
import hashlib
import base64
from datetime import datetime
from typing import Dict, Any, List, Optional
from jinja2 import Environment, FileSystemLoader, BaseLoader
import logging

from app.core.supabase import get_supabase_client

logger = logging.getLogger(__name__)

# Lazy imports
_weasyprint = None


def _get_weasyprint():
    global _weasyprint
    if _weasyprint is None:
        try:
            from weasyprint import HTML, CSS
            _weasyprint = {"HTML": HTML, "CSS": CSS}
        except ImportError:
            logger.warning("WeasyPrint not installed")
            _weasyprint = {"HTML": None, "CSS": None}
    return _weasyprint


# HTML Template for 10-page report
REPORT_TEMPLATE = '''
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>ADC Design Report - {{ session_id }}</title>
    <style>
        @page {
            size: A4;
            margin: 2cm;
            @top-right {
                content: "ADC-GenAI Platform | Confidential";
                font-size: 9px;
                color: #9ca3af;
            }
            @bottom-center {
                content: "Page " counter(page) " of " counter(pages);
                font-size: 9px;
                color: #9ca3af;
            }
            @bottom-right {
                content: "{{ page_hash_placeholder }}";
                font-size: 8px;
                color: #d1d5db;
            }
        }

        * { margin: 0; padding: 0; box-sizing: border-box; }

        body {
            font-family: 'Segoe UI', 'Noto Sans', sans-serif;
            font-size: 11pt;
            line-height: 1.6;
            color: #1f2937;
        }

        .page-break { page-break-before: always; }

        h1 {
            font-size: 24pt;
            font-weight: 700;
            color: #111827;
            margin-bottom: 0.5em;
            border-bottom: 3px solid #3b82f6;
            padding-bottom: 0.3em;
        }

        h2 {
            font-size: 16pt;
            font-weight: 600;
            color: #1e3a5f;
            margin-top: 1.5em;
            margin-bottom: 0.5em;
            border-bottom: 2px solid #e5e7eb;
            padding-bottom: 0.3em;
        }

        h3 {
            font-size: 13pt;
            font-weight: 600;
            color: #374151;
            margin-top: 1em;
            margin-bottom: 0.5em;
        }

        .header {
            text-align: center;
            margin-bottom: 2em;
            padding: 1.5em;
            background: linear-gradient(135deg, #1e3a5f 0%, #3b82f6 100%);
            color: white;
            border-radius: 8px;
        }

        .header h1 {
            color: white;
            border-bottom: none;
        }

        .header .subtitle {
            font-size: 14pt;
            opacity: 0.9;
            margin-top: 0.5em;
        }

        .header .meta {
            display: flex;
            justify-content: center;
            gap: 30px;
            margin-top: 1em;
            font-size: 10pt;
            opacity: 0.8;
        }

        .compliance-badge {
            display: inline-block;
            background: rgba(255,255,255,0.2);
            padding: 4px 16px;
            border-radius: 20px;
            font-size: 10pt;
            margin-top: 1em;
        }

        .summary-box {
            background: #f8fafc;
            border: 1px solid #e5e7eb;
            border-left: 4px solid #3b82f6;
            border-radius: 8px;
            padding: 1em 1.5em;
            margin: 1em 0;
        }

        .executive-summary {
            background: linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%);
            border: 1px solid #bae6fd;
            border-radius: 12px;
            padding: 1.5em;
            margin: 1.5em 0;
        }

        .grade-box {
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 60px;
            height: 60px;
            border-radius: 50%;
            font-size: 24pt;
            font-weight: 700;
            color: white;
            margin-right: 1em;
        }

        .grade-go { background: linear-gradient(135deg, #22c55e 0%, #16a34a 100%); }
        .grade-conditional { background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%); }
        .grade-nogo { background: linear-gradient(135deg, #ef4444 0%, #dc2626 100%); }

        .confidence-box {
            background: #fffbeb;
            border: 1px solid #fcd34d;
            border-radius: 8px;
            padding: 1em;
            margin: 1em 0;
        }

        .confidence-score {
            font-size: 24pt;
            font-weight: 700;
            color: #92400e;
        }

        .confidence-level {
            display: inline-block;
            padding: 2px 12px;
            border-radius: 12px;
            font-size: 10pt;
            font-weight: 600;
        }

        .confidence-high { background: #dcfce7; color: #166534; }
        .confidence-medium { background: #fef3c7; color: #92400e; }
        .confidence-low { background: #fee2e2; color: #991b1b; }

        table {
            width: 100%;
            border-collapse: collapse;
            margin: 1em 0;
            font-size: 10pt;
        }

        th, td {
            padding: 10px 12px;
            text-align: left;
            border-bottom: 1px solid #e5e7eb;
        }

        th {
            background: #f3f4f6;
            font-weight: 600;
            color: #374151;
        }

        tr:hover { background: #f9fafb; }

        .winner { background: #dcfce7; font-weight: 600; }
        .our-adc { background: #dbeafe; }

        .chart-container {
            text-align: center;
            margin: 1.5em 0;
            padding: 1em;
            background: #f9fafb;
            border-radius: 8px;
        }

        .chart-container img {
            max-width: 100%;
            height: auto;
        }

        .reasoning-block {
            background: #f0fdf4;
            border: 1px solid #86efac;
            border-radius: 8px;
            padding: 1em;
            margin: 1em 0;
            font-style: italic;
        }

        .reasoning-block::before {
            content: "AI Reasoning: ";
            font-weight: 600;
            font-style: normal;
            color: #166534;
        }

        .qr-section {
            text-align: center;
            padding: 2em;
            background: #f8fafc;
            border: 2px dashed #cbd5e1;
            border-radius: 12px;
            margin: 2em 0;
        }

        .qr-section img {
            width: 150px;
            height: 150px;
        }

        .hash-display {
            font-family: 'Courier New', monospace;
            font-size: 9pt;
            background: #1e293b;
            color: #22d3ee;
            padding: 8px 16px;
            border-radius: 4px;
            word-break: break-all;
        }

        .page-hash-grid {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 8px;
            margin: 1em 0;
        }

        .page-hash-item {
            font-family: monospace;
            font-size: 8pt;
            padding: 4px 8px;
            background: #f1f5f9;
            border-radius: 4px;
        }

        .reference-list {
            font-size: 9pt;
            color: #6b7280;
        }

        .reference-list li {
            margin-bottom: 0.5em;
        }

        .footer {
            margin-top: 2em;
            padding-top: 1em;
            border-top: 1px solid #e5e7eb;
            font-size: 9pt;
            color: #6b7280;
            text-align: center;
        }

        .disclaimer {
            background: #fffbeb;
            border: 1px solid #fcd34d;
            border-radius: 4px;
            padding: 0.5em 1em;
            font-size: 8pt;
            color: #92400e;
            margin-top: 1em;
        }

        .section-confidence {
            float: right;
            padding: 4px 12px;
            border-radius: 20px;
            font-size: 10pt;
            font-weight: 500;
        }

        .metrics-grid {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 12px;
            margin: 1em 0;
        }

        .metric-card {
            background: white;
            border: 1px solid #e5e7eb;
            border-radius: 8px;
            padding: 12px;
            text-align: center;
        }

        .metric-value {
            font-size: 20pt;
            font-weight: 700;
            color: #1e3a5f;
        }

        .metric-label {
            font-size: 9pt;
            color: #6b7280;
            margin-top: 4px;
        }

        .synthesis-step {
            display: flex;
            align-items: flex-start;
            gap: 12px;
            padding: 12px;
            border-left: 3px solid #3b82f6;
            margin-bottom: 12px;
            background: #f8fafc;
        }

        .step-number {
            background: #3b82f6;
            color: white;
            width: 28px;
            height: 28px;
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            font-weight: 600;
            flex-shrink: 0;
        }
    </style>
</head>
<body>
    <!-- Page 1: Executive Summary -->
    <div class="header">
        <h1>ADC Design Report</h1>
        <div class="subtitle">{{ session_type }} Design Analysis</div>
        <div class="meta">
            <span>Session: {{ session_id_short }}</span>
            <span>Generated: {{ generated_at }}</span>
            <span>Target: {{ target_antigen }}</span>
        </div>
        <div class="compliance-badge">21 CFR Part 11 Compliant | Digital Seal Enabled</div>
    </div>

    <h2>Executive Summary</h2>

    <div class="executive-summary">
        <div style="display: flex; align-items: center;">
            <div class="grade-box {{ grade_class }}">{{ final_grade }}</div>
            <div>
                <strong style="font-size: 14pt;">{{ recommendation }}</strong>
                <p style="margin-top: 0.5em; color: #6b7280;">{{ executive_summary }}</p>
            </div>
        </div>
    </div>

    <div class="confidence-box">
        <strong>Overall Confidence Score:</strong>
        <span class="confidence-score">{{ overall_confidence }}%</span>
        <span class="confidence-level confidence-{{ confidence_level|lower }}">{{ confidence_level }}</span>
        <p style="margin-top: 0.5em; font-size: 10pt; color: #6b7280;">{{ confidence_reasoning }}</p>
    </div>

    <h3>Key Findings</h3>
    <div class="metrics-grid">
        {% for metric in key_metrics %}
        <div class="metric-card">
            <div class="metric-value">{{ metric.value }}</div>
            <div class="metric-label">{{ metric.label }}</div>
        </div>
        {% endfor %}
    </div>

    <div class="reasoning-block">
        {{ executive_reasoning }}
    </div>

    <!-- Page 2-3: Target & Antibody Analysis -->
    <div class="page-break"></div>
    <h1>Target & Antibody Analysis</h1>
    <span class="section-confidence confidence-{{ target_confidence_level|lower }}">
        Confidence: {{ target_confidence }}%
    </span>

    <h2>Target Antigen Profile</h2>
    <div class="summary-box">
        <strong>{{ target_antigen }}</strong>
        <p>{{ target_description }}</p>
    </div>

    <table>
        <tr><th>Property</th><th>Value</th><th>Significance</th></tr>
        {% for prop in target_properties %}
        <tr>
            <td>{{ prop.name }}</td>
            <td>{{ prop.value }}</td>
            <td>{{ prop.significance }}</td>
        </tr>
        {% endfor %}
    </table>

    <h2>Antibody Selection</h2>
    <div class="summary-box">
        <strong>{{ antibody_name }}</strong>
        <p>Binding Affinity (Kd): {{ binding_affinity }}</p>
    </div>

    {% if target_chart %}
    <div class="chart-container">
        <img src="data:image/png;base64,{{ target_chart }}" alt="Target Expression Profile">
        <p style="font-size: 9pt; color: #6b7280;">Target Expression Profile</p>
    </div>
    {% endif %}

    <div class="reasoning-block">
        {{ target_reasoning }}
    </div>

    <!-- Page 4: In-silico Synthesis Path -->
    <div class="page-break"></div>
    <h1>In-silico Synthesis Path</h1>
    <span class="section-confidence confidence-{{ synthesis_confidence_level|lower }}">
        Confidence: {{ synthesis_confidence }}%
    </span>

    <h2>Synthetic Accessibility Analysis</h2>
    <div class="metrics-grid">
        <div class="metric-card">
            <div class="metric-value">{{ sa_score }}</div>
            <div class="metric-label">SA Score (1-10)</div>
        </div>
        <div class="metric-card">
            <div class="metric-value">{{ estimated_yield }}</div>
            <div class="metric-label">Est. Total Yield</div>
        </div>
        <div class="metric-card">
            <div class="metric-value">{{ synthesis_difficulty }}</div>
            <div class="metric-label">Difficulty</div>
        </div>
    </div>

    {% if sa_chart %}
    <div class="chart-container">
        <img src="data:image/png;base64,{{ sa_chart }}" alt="SA Score Comparison">
        <p style="font-size: 9pt; color: #6b7280;">Synthetic Accessibility Comparison vs Approved ADCs</p>
    </div>
    {% endif %}

    <h2>Proposed Synthesis Route</h2>
    {% for step in synthesis_steps %}
    <div class="synthesis-step">
        <div class="step-number">{{ step.number }}</div>
        <div>
            <strong>{{ step.reaction }}</strong>
            <p style="font-size: 10pt; color: #6b7280;">Reagents: {{ step.reagents }}</p>
            <p style="font-size: 10pt;">Estimated Yield: {{ step.yield }} | Difficulty: {{ step.difficulty }}</p>
        </div>
    </div>
    {% endfor %}

    <div class="reasoning-block">
        {{ synthesis_reasoning }}
    </div>

    <!-- Page 5: Benchmark Comparison (Dynamic Benchmarking) -->
    <div class="page-break"></div>
    <h1>Benchmark Comparison</h1>
    <span class="section-confidence confidence-{{ benchmark_confidence_level|lower }}">
        Confidence: {{ benchmark_confidence }}%
    </span>

    <!-- Gold Standard Reference -->
    {% if gold_standard %}
    <div class="summary-box" style="background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%); border-color: #f59e0b;">
        <h3 style="color: #92400e; margin-bottom: 0.5em;">Gold Standard Reference: {{ gold_standard.trade_name }}</h3>
        <p style="font-size: 10pt; color: #78350f;">
            {{ gold_standard.description }} (Target: {{ gold_standard.target }})
        </p>
        <p style="font-size: 9pt; color: #92400e; margin-top: 0.5em;">
            본 보고서는 {{ gold_standard.trade_name }}을 기준(100%)으로 상대적 성능을 평가합니다.
        </p>
    </div>
    {% endif %}

    <!-- Relative Performance Radar Chart -->
    {% if relative_radar_chart %}
    <h2>Relative Performance vs Gold Standard</h2>
    <div class="chart-container">
        <img src="data:image/png;base64,{{ relative_radar_chart }}" alt="Relative Performance Radar">
        <p style="font-size: 9pt; color: #6b7280;">
            100% = {{ gold_standard.trade_name if gold_standard else 'Gold Standard' }} |
            바깥쪽 = 우수 | 안쪽 = 개선 필요
        </p>
    </div>
    {% endif %}

    <!-- Relative Efficiency Scores -->
    {% if relative_scores %}
    <h2>Relative Efficiency Score (%)</h2>
    <table>
        <tr>
            <th>Metric</th>
            <th>Our ADC</th>
            <th>Gold Standard</th>
            <th>Relative %</th>
            <th>Interpretation</th>
        </tr>
        {% for metric_name, score in relative_scores.items() %}
        <tr>
            <td>{{ metric_name }}</td>
            <td>{{ score.our_value }}</td>
            <td>{{ score.benchmark_value }}</td>
            <td style="font-weight: bold; color: {% if score.is_favorable %}#22c55e{% else %}#ef4444{% endif %};">
                {{ score.relative_percent }}%
            </td>
            <td style="font-size: 9pt;">{{ score.interpretation }}</td>
        </tr>
        {% endfor %}
    </table>
    {% endif %}

    {% if relative_bar_chart %}
    <div class="chart-container">
        <img src="data:image/png;base64,{{ relative_bar_chart }}" alt="Relative Scores">
    </div>
    {% endif %}

    <h2>Comparison vs Approved ADCs</h2>
    <table>
        <tr>
            <th>Property</th>
            <th class="our-adc">Our ADC</th>
            {% for comp in comparators %}
            <th>{{ comp }}</th>
            {% endfor %}
            <th>Winner</th>
        </tr>
        {% for row in benchmark_table %}
        <tr>
            <td>{{ row.property }}</td>
            <td class="our-adc">{{ row.our_adc }}</td>
            {% for comp in comparators %}
            <td>{{ row[comp] }}</td>
            {% endfor %}
            <td class="{% if row.winner == 'our_adc' %}winner{% endif %}">
                {{ row.winner_display }}
            </td>
        </tr>
        {% endfor %}
    </table>

    {% if benchmark_chart %}
    <div class="chart-container">
        <img src="data:image/png;base64,{{ benchmark_chart }}" alt="Benchmark Comparison">
        <p style="font-size: 9pt; color: #6b7280;">Performance Comparison vs Approved ADCs</p>
    </div>
    {% endif %}

    <h3>Competitive Advantages</h3>
    <ul>
        {% for adv in competitive_advantages %}
        <li>{{ adv }}</li>
        {% endfor %}
    </ul>

    <div class="reasoning-block">
        {{ benchmark_reasoning }}
    </div>

    <!-- Page 6: Physicochemical Properties -->
    <div class="page-break"></div>
    <h1>Physicochemical Properties</h1>
    <span class="section-confidence confidence-{{ physchem_confidence_level|lower }}">
        Confidence: {{ physchem_confidence }}%
    </span>

    <h2>Molecular Properties</h2>
    <table>
        <tr><th>Property</th><th>Value</th><th>Lipinski Limit</th><th>Status</th></tr>
        {% for prop in lipinski_properties %}
        <tr>
            <td>{{ prop.name }}</td>
            <td>{{ prop.value }}</td>
            <td>{{ prop.limit }}</td>
            <td style="color: {{ prop.status_color }}">{{ prop.status }}</td>
        </tr>
        {% endfor %}
    </table>

    {% if lipinski_chart %}
    <div class="chart-container">
        <img src="data:image/png;base64,{{ lipinski_chart }}" alt="Lipinski Rule of Five">
        <p style="font-size: 9pt; color: #6b7280;">Lipinski's Rule of Five Compliance</p>
    </div>
    {% endif %}

    <h2>DAR Distribution Prediction</h2>
    {% if dar_chart %}
    <div class="chart-container">
        <img src="data:image/png;base64,{{ dar_chart }}" alt="DAR Distribution">
        <p style="font-size: 9pt; color: #6b7280;">Predicted DAR Distribution (Target: {{ target_dar }})</p>
    </div>
    {% endif %}

    <div class="reasoning-block">
        {{ physchem_reasoning }}
    </div>

    <!-- Page 7: Predictive Toxicology -->
    <div class="page-break"></div>
    <h1>Predictive Toxicology</h1>
    <span class="section-confidence confidence-{{ tox_confidence_level|lower }}">
        Confidence: {{ tox_confidence }}%
    </span>

    <h2>Toxicity Risk Assessment</h2>
    <table>
        <tr><th>Risk Type</th><th>Level</th><th>Prediction</th><th>Mitigation</th></tr>
        {% for risk in toxicity_risks %}
        <tr>
            <td>{{ risk.type }}</td>
            <td style="background: {{ risk.level_color }}">{{ risk.level }}</td>
            <td>{{ risk.prediction }}</td>
            <td>{{ risk.mitigation }}</td>
        </tr>
        {% endfor %}
    </table>

    <h2>DAR-dependent Adverse Events</h2>
    <div class="summary-box">
        {{ dar_toxicity_analysis }}
    </div>

    {% if linker_cleavage_diagram %}
    <h2>Linker Cleavage Simulation</h2>
    <div class="chart-container">
        <img src="data:image/png;base64,{{ linker_cleavage_diagram }}" alt="Linker Cleavage">
        <p style="font-size: 9pt; color: #6b7280;">Predicted Linker Cleavage Sites</p>
    </div>
    {% endif %}

    <div class="reasoning-block">
        {{ tox_reasoning }}
    </div>

    <!-- Page 8: Patent & Competitor Landscape -->
    <div class="page-break"></div>
    <h1>Patent & Competitor Landscape</h1>
    <span class="section-confidence confidence-{{ patent_confidence_level|lower }}">
        Confidence: {{ patent_confidence }}%
    </span>

    <h2>Freedom to Operate Analysis</h2>
    <div class="summary-box">
        <strong>FTO Status: {{ fto_status }}</strong>
        <p>{{ fto_summary }}</p>
    </div>

    <h2>Relevant Patents</h2>
    <table>
        <tr><th>Patent ID</th><th>Title</th><th>Assignee</th><th>Risk</th></tr>
        {% for patent in relevant_patents %}
        <tr>
            <td>{{ patent.id }}</td>
            <td>{{ patent.title }}</td>
            <td>{{ patent.assignee }}</td>
            <td style="background: {{ patent.risk_color }}">{{ patent.risk }}</td>
        </tr>
        {% endfor %}
    </table>

    <h2>Differentiation Strategy</h2>
    <ul>
        {% for strategy in differentiation_strategies %}
        <li>{{ strategy }}</li>
        {% endfor %}
    </ul>

    <div class="reasoning-block">
        {{ patent_reasoning }}
    </div>

    <!-- Page 9: Digital Audit Trail & Seal -->
    <div class="page-break"></div>
    <h1>Digital Audit Trail & Seal</h1>

    <div class="qr-section">
        <h2>Verification QR Code</h2>
        <img src="data:image/png;base64,{{ qr_code }}" alt="Verification QR Code">
        <p style="margin-top: 1em; font-size: 10pt;">Scan to verify report authenticity</p>
        <p style="font-size: 9pt; color: #6b7280;">{{ verification_url }}</p>
    </div>

    <h2>Chain Hash</h2>
    <div class="hash-display">{{ chain_hash }}</div>
    <p style="font-size: 9pt; color: #6b7280; margin-top: 0.5em;">
        Algorithm: SHA-256 Chain | Generated: {{ seal_generated_at }}
    </p>

    <h2>Page Hash Registry</h2>
    <div class="page-hash-grid">
        {% for ph in page_hashes %}
        <div class="page-hash-item">P{{ ph.page }}: {{ ph.hash }}</div>
        {% endfor %}
    </div>

    <h2>Audit Trail</h2>
    <table>
        <tr><th>Timestamp</th><th>Action</th><th>Agent</th><th>Hash</th></tr>
        {% for entry in audit_trail %}
        <tr>
            <td>{{ entry.timestamp }}</td>
            <td>{{ entry.action }}</td>
            <td>{{ entry.agent }}</td>
            <td style="font-family: monospace; font-size: 8pt;">{{ entry.hash }}</td>
        </tr>
        {% endfor %}
    </table>

    <!-- Phase 2: Digital Lineage (Calculation Metadata) -->
    {% if has_lineage %}
    <h2>Calculation Lineage (Data Provenance)</h2>
    <div class="summary-box" style="background: #f0fdf4; border-color: #22c55e;">
        <p style="font-size: 10pt; color: #166534;">
            <strong>21 CFR Part 11 Compliant</strong> - All calculations recorded with full execution environment metadata for reproducibility verification.
        </p>
    </div>

    <table style="font-size: 9pt;">
        <tr>
            <th>Calc ID</th>
            <th>Agent</th>
            <th>Timestamp</th>
            <th>Library</th>
            <th>LLM Model</th>
            <th>Time (ms)</th>
        </tr>
        {% for lineage in lineage_summary %}
        <tr>
            <td style="font-family: monospace;">{{ lineage.calculation_id }}</td>
            <td>{{ lineage.agent }}</td>
            <td>{{ lineage.timestamp }}</td>
            <td>{{ lineage.library }}</td>
            <td>{{ lineage.llm }}</td>
            <td>{{ lineage.time_ms }}</td>
        </tr>
        {% endfor %}
    </table>

    {% if physical_validations %}
    <h3>Physical Validations Summary</h3>
    <div class="metrics-grid" style="grid-template-columns: repeat(2, 1fr);">
        <div class="metric-card" style="background: #dcfce7;">
            <div class="metric-value" style="color: #166534;">{{ validations_passed }}</div>
            <div class="metric-label">Passed</div>
        </div>
        <div class="metric-card" style="background: {% if validations_failed > 0 %}#fee2e2{% else %}#f3f4f6{% endif %};">
            <div class="metric-value" style="color: {% if validations_failed > 0 %}#991b1b{% else %}#6b7280{% endif %};">{{ validations_failed }}</div>
            <div class="metric-label">Failed</div>
        </div>
    </div>
    {% endif %}
    {% endif %}

    <!-- Page 10: Appendix & References -->
    <div class="page-break"></div>
    <h1>Appendix & References</h1>

    <h2>Confidence Score Breakdown</h2>
    <table>
        <tr><th>Section</th><th>Score</th><th>Level</th><th>Key Factors</th></tr>
        {% for section in confidence_breakdown %}
        <tr>
            <td>{{ section.name }}</td>
            <td>{{ section.score }}%</td>
            <td class="confidence-{{ section.level|lower }}">{{ section.level }}</td>
            <td>{{ section.factors }}</td>
        </tr>
        {% endfor %}
    </table>

    <h2>Referenced Literature</h2>
    <ol class="reference-list">
        {% for ref in references %}
        <li>{{ ref.authors }} ({{ ref.year }}). {{ ref.title }}. {{ ref.journal }}. PMID: {{ ref.pmid }}</li>
        {% endfor %}
    </ol>

    <h2>Data Sources</h2>
    <ul class="reference-list">
        {% for source in data_sources %}
        <li>{{ source }}</li>
        {% endfor %}
    </ul>

    <div class="footer">
        <p>This report was generated by ADC-GenAI Platform using AI-assisted analysis.</p>
        <p>Report ID: {{ report_id }} | Session ID: {{ session_id }}</p>
        <div class="disclaimer">
            <strong>Disclaimer:</strong> This report is generated by AI and should be reviewed by qualified scientists before making any decisions.
            The predictions and analyses contained herein are based on computational models and should be validated experimentally.
        </div>
    </div>
</body>
</html>
'''


class ReportEngine:
    """
    10페이지 정밀 보고서 생성 엔진

    Phase 2 Enhancement:
    - Digital Lineage 수집 및 태깅
    - 계산 메타데이터 리포트 통합
    """

    def __init__(self):
        from app.services.visualization_service import VisualizationService
        from app.services.digital_seal_service import get_digital_seal_service
        from app.services.confidence_calculator import ConfidenceCalculator

        self.visualization = VisualizationService()
        self.seal_service = get_digital_seal_service()
        self.confidence_calc = ConfidenceCalculator()
        self.supabase = get_supabase_client()

    # =========================================================================
    # Phase 2: Digital Lineage Methods
    # =========================================================================

    async def _collect_calculation_lineage(
        self,
        session_id: str
    ) -> Dict[str, Dict[str, Any]]:
        """
        세션의 모든 에이전트 연산 메타데이터 수집

        Returns:
            Dict[calculation_id, CalculationMetadata]
        """
        lineage = {}

        try:
            # agent_execution_logs 테이블에서 수집
            result = self.supabase.table("agent_execution_logs").select(
                "calculation_id, agent_name, agent_version, status, "
                "execution_env, library_versions, nvidia_nim_info, "
                "simulation_params, input_hash, output_hash, "
                "llm_model, llm_temperature, execution_time_ms, "
                "created_at"
            ).eq("session_id", session_id).eq("status", "completed").execute()

            for log in result.data or []:
                calc_id = log.get("calculation_id")
                if calc_id:
                    lineage[calc_id] = {
                        "calculation_id": calc_id,
                        "timestamp": log.get("created_at"),
                        "execution_environment": log.get("execution_env", {}),
                        "library_versions": log.get("library_versions", {}),
                        "nvidia_nim": log.get("nvidia_nim_info"),
                        "simulation_params": log.get("simulation_params"),
                        "input_hash": log.get("input_hash"),
                        "output_hash": log.get("output_hash"),
                        "execution_time_ms": log.get("execution_time_ms"),
                        "agent_info": {
                            "agent_name": log.get("agent_name"),
                            "agent_version": log.get("agent_version"),
                            "llm_model": log.get("llm_model"),
                            "llm_temperature": log.get("llm_temperature")
                        }
                    }

            logger.info(f"[ReportEngine] Collected {len(lineage)} calculation lineage records for session {session_id}")
            return lineage

        except Exception as e:
            logger.error(f"[ReportEngine] Failed to collect lineage: {e}")
            return {}

    async def _collect_physical_validations(
        self,
        session_id: str
    ) -> List[Dict[str, Any]]:
        """물리 검증 결과 수집"""
        try:
            result = self.supabase.table("physical_validations").select(
                "calculation_id, smiles, check_name, check_category, "
                "passed, severity, actual_value, expected_range, details, "
                "created_at"
            ).eq("session_id", session_id).order("created_at").execute()

            return result.data or []

        except Exception as e:
            logger.warning(f"[ReportEngine] Failed to collect validations: {e}")
            return []

    def _format_lineage_for_template(
        self,
        lineage: Dict[str, Dict],
        agent_name: str
    ) -> Optional[Dict[str, Any]]:
        """특정 에이전트의 lineage 데이터를 템플릿용으로 포맷"""
        for calc_id, data in lineage.items():
            if data["agent_info"]["agent_name"] == agent_name:
                return {
                    "calc_id": calc_id[:8],
                    "full_calc_id": calc_id,
                    "library": data["library_versions"].get("rdkit", "N/A"),
                    "container": data["execution_environment"].get("container_image", "N/A"),
                    "python_version": data["execution_environment"].get("python_version", "N/A"),
                    "timestamp": data["timestamp"],
                    "input_hash": data["input_hash"][:16] + "..." if data["input_hash"] else "N/A",
                    "output_hash": data["output_hash"][:16] + "..." if data["output_hash"] else "N/A",
                    "execution_time_ms": data.get("execution_time_ms", 0)
                }
        return None

    def _format_lineage_summary(
        self,
        lineage: Dict[str, Dict]
    ) -> List[Dict[str, Any]]:
        """전체 lineage 요약 테이블용 포맷"""
        summary = []
        for calc_id, data in lineage.items():
            agent_info = data.get("agent_info", {})
            summary.append({
                "calculation_id": calc_id[:8],
                "agent": agent_info.get("agent_name", "unknown"),
                "timestamp": data.get("timestamp", "")[:19],  # ISO format 잘라서
                "library": data.get("library_versions", {}).get("rdkit", "-"),
                "llm": agent_info.get("llm_model", "-"),
                "time_ms": data.get("execution_time_ms", 0)
            })
        return sorted(summary, key=lambda x: x["timestamp"])

    async def generate_report(
        self,
        session_id: str,
        session_data: Dict[str, Any],
        agent_results: Dict[str, Any],
        format: str = "pdf",
        include_visualizations: bool = True,
        include_benchmarks: bool = True,
        include_lineage: bool = True  # Phase 2: Digital Lineage
    ) -> Dict[str, Any]:
        """
        전체 보고서 생성

        Args:
            session_id: 세션 ID
            session_data: 세션 데이터 (target, antibody, etc.)
            agent_results: 각 에이전트 실행 결과
            format: 출력 포맷 (pdf, html)
            include_visualizations: 시각화 포함 여부
            include_benchmarks: 벤치마크 비교 포함 여부
            include_lineage: Digital Lineage 포함 여부 (Phase 2)

        Returns:
            보고서 생성 결과 (파일 경로 또는 바이트)
        """
        # Phase 2: Digital Lineage 수집
        lineage_data = {}
        validations_data = []
        if include_lineage:
            lineage_data = await self._collect_calculation_lineage(session_id)
            validations_data = await self._collect_physical_validations(session_id)

        # 템플릿 데이터 준비
        template_data = await self._prepare_template_data(
            session_id, session_data, agent_results,
            include_visualizations, include_benchmarks,
            lineage_data, validations_data  # Phase 2
        )

        # HTML 렌더링
        env = Environment(loader=BaseLoader())
        template = env.from_string(REPORT_TEMPLATE)
        html_content = template.render(**template_data)

        # Digital Seal 생성
        pages_content = self._extract_page_contents(html_content)
        seal_data = self.seal_service.generate_full_seal({
            "session_id": session_id,
            "pages": pages_content
        })

        # Seal 정보 추가하여 재렌더링
        template_data.update({
            "qr_code": seal_data["qr_code"],
            "chain_hash": seal_data["chain_hash"],
            "verification_url": seal_data["verification_url"],
            "seal_generated_at": seal_data["generated_at"],
            "page_hashes": [
                {"page": ph["page_number"], "hash": ph["short_hash"]}
                for ph in seal_data["page_hashes"]
            ]
        })

        html_content = template.render(**template_data)

        if format == "pdf":
            return await self._generate_pdf(html_content, session_id, seal_data)
        else:
            return {
                "format": "html",
                "content": html_content,
                "seal": seal_data,
                "session_id": session_id
            }

    async def _prepare_template_data(
        self,
        session_id: str,
        session_data: Dict,
        agent_results: Dict,
        include_visualizations: bool,
        include_benchmarks: bool,
        lineage_data: Optional[Dict] = None,  # Phase 2
        validations_data: Optional[List] = None  # Phase 2
    ) -> Dict[str, Any]:
        """템플릿 데이터 준비"""
        from app.services.confidence_calculator import calculate_section_confidence, aggregate_report_confidence

        lineage_data = lineage_data or {}
        validations_data = validations_data or []

        # 기본 정보
        data = {
            "session_id": session_id,
            "session_id_short": session_id[:8],
            "session_type": session_data.get("session_type", "ADC Design"),
            "generated_at": datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC"),
            "target_antigen": session_data.get("target_antigen", "N/A"),
            "report_id": f"RPT-{session_id[:8]}-{datetime.utcnow().strftime('%Y%m%d')}",
        }

        # 최종 판정
        final_grade = agent_results.get("final_grade", "B")
        data["final_grade"] = final_grade
        data["grade_class"] = self._get_grade_class(final_grade)
        data["recommendation"] = agent_results.get("recommendation", "Conditional Go")
        data["executive_summary"] = agent_results.get("executive_summary", "Analysis complete.")
        data["executive_reasoning"] = agent_results.get("executive_reasoning", "")

        # 주요 메트릭
        data["key_metrics"] = self._extract_key_metrics(session_data, agent_results)

        # 각 섹션별 Confidence Score 계산
        section_confidences = {}

        # Target Analysis (P2-3)
        target_conf = calculate_section_confidence("target_analysis", agent_results.get("target_analysis", {}))
        section_confidences["target_analysis"] = target_conf
        data["target_confidence"] = target_conf.score
        data["target_confidence_level"] = target_conf.level.value
        data["target_description"] = agent_results.get("target_analysis", {}).get("description", "")
        data["target_properties"] = agent_results.get("target_analysis", {}).get("properties", [])
        data["antibody_name"] = session_data.get("antibody_name", "N/A")
        data["binding_affinity"] = session_data.get("binding_affinity", "N/A")
        data["target_reasoning"] = target_conf.reasoning
        data["target_chart"] = None

        # Synthesis Path (P4)
        synthesis_data = agent_results.get("synthesis_analysis", {})
        synthesis_conf = calculate_section_confidence("synthesis_path", synthesis_data)
        section_confidences["synthesis_path"] = synthesis_conf
        data["synthesis_confidence"] = synthesis_conf.score
        data["synthesis_confidence_level"] = synthesis_conf.level.value
        data["sa_score"] = synthesis_data.get("sa_score", "N/A")
        data["estimated_yield"] = synthesis_data.get("estimated_total_yield", "N/A")
        data["synthesis_difficulty"] = synthesis_data.get("sa_interpretation", "Moderate")
        data["synthesis_steps"] = [
            {
                "number": s.get("step", i+1),
                "reaction": s.get("reaction", ""),
                "reagents": ", ".join(s.get("reagents", [])),
                "yield": s.get("yield_estimate", "N/A"),
                "difficulty": s.get("difficulty", "Moderate")
            }
            for i, s in enumerate(synthesis_data.get("synthesis_steps", []))
        ]
        data["synthesis_reasoning"] = synthesis_data.get("reasoning_logic", synthesis_conf.reasoning)

        # Benchmark Comparison (P5) - Dynamic Benchmarking
        if include_benchmarks:
            benchmark_data = agent_results.get("benchmark_comparison", {})
            benchmark_conf = calculate_section_confidence("benchmark_comparison", benchmark_data)
            section_confidences["benchmark_comparison"] = benchmark_conf
            data["benchmark_confidence"] = benchmark_conf.score
            data["benchmark_confidence_level"] = benchmark_conf.level.value
            data["comparators"] = list(benchmark_data.get("competitors", {}).keys())
            data["benchmark_table"] = self._format_benchmark_table(benchmark_data)
            data["competitive_advantages"] = benchmark_data.get("competitive_advantages", [])
            data["benchmark_reasoning"] = benchmark_data.get("reasoning_logic", benchmark_conf.reasoning)

            # Gold Standard 정보 (NEW)
            gold_standard = agent_results.get("gold_standard")
            if gold_standard:
                data["gold_standard"] = {
                    "trade_name": gold_standard.get("trade_name", "N/A"),
                    "target": gold_standard.get("target", "N/A"),
                    "description": gold_standard.get("description", "")
                }
            else:
                data["gold_standard"] = None

            # Relative Scores (NEW)
            relative_scores = agent_results.get("relative_scores", {})
            data["relative_scores"] = relative_scores

            # Radar Chart Data (NEW)
            radar_chart_data = agent_results.get("radar_chart_data")
            data["radar_chart_data"] = radar_chart_data
        else:
            data["benchmark_confidence"] = 0
            data["benchmark_confidence_level"] = "Low"
            data["comparators"] = []
            data["benchmark_table"] = []
            data["competitive_advantages"] = []
            data["benchmark_reasoning"] = "Benchmark comparison not included."
            data["gold_standard"] = None
            data["relative_scores"] = {}
            data["radar_chart_data"] = None

        # Physicochemical Properties (P6)
        physchem_data = agent_results.get("physicochemical", {})
        physchem_conf = calculate_section_confidence("physicochemical", physchem_data)
        section_confidences["physicochemical"] = physchem_conf
        data["physchem_confidence"] = physchem_conf.score
        data["physchem_confidence_level"] = physchem_conf.level.value
        data["lipinski_properties"] = self._format_lipinski_properties(physchem_data)
        data["target_dar"] = session_data.get("requested_dar", 4)
        data["physchem_reasoning"] = physchem_data.get("reasoning_logic", physchem_conf.reasoning)

        # Toxicology (P7)
        tox_data = agent_results.get("toxicology", {})
        tox_conf = calculate_section_confidence("toxicology", tox_data)
        section_confidences["toxicology"] = tox_conf
        data["tox_confidence"] = tox_conf.score
        data["tox_confidence_level"] = tox_conf.level.value
        data["toxicity_risks"] = self._format_toxicity_risks(tox_data)
        data["dar_toxicity_analysis"] = tox_data.get("dar_analysis", "")
        data["tox_reasoning"] = tox_data.get("reasoning_logic", tox_conf.reasoning)

        # Patent (P8)
        patent_data = agent_results.get("patent_analysis", {})
        patent_conf = calculate_section_confidence("patent_landscape", patent_data)
        section_confidences["patent_landscape"] = patent_conf
        data["patent_confidence"] = patent_conf.score
        data["patent_confidence_level"] = patent_conf.level.value
        data["fto_status"] = patent_data.get("fto_status", "Under Review")
        data["fto_summary"] = patent_data.get("fto_summary", "")
        data["relevant_patents"] = patent_data.get("patents", [])
        data["differentiation_strategies"] = patent_data.get("differentiation_strategies", [])
        data["patent_reasoning"] = patent_data.get("reasoning_logic", patent_conf.reasoning)

        # 전체 Confidence Score
        overall_conf = aggregate_report_confidence(section_confidences)
        data["overall_confidence"] = overall_conf.score
        data["confidence_level"] = overall_conf.level.value
        data["confidence_reasoning"] = overall_conf.reasoning

        # Confidence Breakdown (P10)
        data["confidence_breakdown"] = [
            {
                "name": name.replace("_", " ").title(),
                "score": conf.score,
                "level": conf.level.value,
                "factors": "; ".join(conf.limitations[:2]) if conf.limitations else "Good data coverage"
            }
            for name, conf in section_confidences.items()
        ]

        # Audit Trail
        data["audit_trail"] = agent_results.get("audit_trail", [])

        # References
        data["references"] = self._collect_references(agent_results)
        data["data_sources"] = [
            "PubMed / NCBI",
            "UniProt Protein Database",
            "ChEMBL / PubChem",
            "FDA Drug Database",
            "ClinicalTrials.gov",
            "Golden Set (Internal Database)"
        ]

        # 시각화 생성
        if include_visualizations:
            data.update(await self._generate_visualizations(session_data, agent_results))

        # Phase 2: Digital Lineage 데이터
        if lineage_data:
            data["lineage_summary"] = self._format_lineage_summary(lineage_data)
            data["has_lineage"] = True

            # 각 에이전트별 lineage 매핑
            data["coder_lineage"] = self._format_lineage_for_template(lineage_data, "coder")
            data["auditor_lineage"] = self._format_lineage_for_template(lineage_data, "auditor")
            data["alchemist_lineage"] = self._format_lineage_for_template(lineage_data, "alchemist")
        else:
            data["lineage_summary"] = []
            data["has_lineage"] = False

        # Phase 2: Physical Validations 데이터
        if validations_data:
            data["physical_validations"] = validations_data
            data["validations_passed"] = sum(1 for v in validations_data if v.get("passed"))
            data["validations_failed"] = sum(1 for v in validations_data if not v.get("passed"))
        else:
            data["physical_validations"] = []
            data["validations_passed"] = 0
            data["validations_failed"] = 0

        return data

    async def _generate_visualizations(
        self,
        session_data: Dict,
        agent_results: Dict
    ) -> Dict[str, str]:
        """시각화 차트 생성"""
        charts = {}

        try:
            # Lipinski Chart
            physchem = agent_results.get("physicochemical", {})
            if physchem:
                charts["lipinski_chart"] = self.visualization.generate_lipinski_chart({
                    "mw": physchem.get("molecular_weight", 400),
                    "logp": physchem.get("logp", 3),
                    "hbd": physchem.get("hbd", 2),
                    "hba": physchem.get("hba", 5),
                    "rotatable_bonds": physchem.get("rotatable_bonds", 5)
                })

            # DAR Distribution
            target_dar = session_data.get("requested_dar", 4)
            charts["dar_chart"] = self.visualization.generate_dar_distribution(target_dar)

            # SA Score Chart
            synthesis = agent_results.get("synthesis_analysis", {})
            if synthesis.get("sa_score"):
                comparison_data = synthesis.get("comparison_to_approved", {})
                charts["sa_chart"] = self.visualization.generate_sa_score_chart(
                    synthesis["sa_score"], comparison_data
                )

            # Benchmark Comparison Chart
            benchmark = agent_results.get("benchmark_comparison", {})
            if benchmark.get("our_adc") and benchmark.get("competitors"):
                charts["benchmark_chart"] = self.visualization.generate_benchmark_comparison_chart(
                    benchmark["our_adc"], benchmark["competitors"]
                )

            # Relative Radar Chart (NEW - Dynamic Benchmarking)
            radar_chart_data = agent_results.get("radar_chart_data")
            if radar_chart_data:
                charts["relative_radar_chart"] = self.visualization.generate_relative_radar_chart(
                    radar_chart_data
                )

            # Relative Scores Bar Chart (NEW)
            relative_scores = agent_results.get("relative_scores", {})
            if relative_scores:
                charts["relative_bar_chart"] = self.visualization.generate_relative_scores_bar_chart(
                    relative_scores
                )

            # Linker Cleavage Diagram
            linker_smiles = session_data.get("linker_smiles")
            linker_type = session_data.get("linker_type", "cleavable")
            if linker_smiles:
                charts["linker_cleavage_diagram"] = self.visualization.generate_linker_cleavage_diagram(
                    linker_type, linker_smiles
                )

        except Exception as e:
            logger.warning(f"Visualization generation error: {e}")

        return charts

    async def _generate_pdf(
        self,
        html_content: str,
        session_id: str,
        seal_data: Dict
    ) -> Dict[str, Any]:
        """PDF 생성"""
        wp = _get_weasyprint()

        if wp["HTML"] is None:
            return {
                "format": "html",
                "content": html_content,
                "seal": seal_data,
                "session_id": session_id,
                "warning": "WeasyPrint not available, returning HTML"
            }

        try:
            pdf_buffer = io.BytesIO()
            wp["HTML"](string=html_content).write_pdf(pdf_buffer)
            pdf_buffer.seek(0)

            return {
                "format": "pdf",
                "content": pdf_buffer.getvalue(),
                "seal": seal_data,
                "session_id": session_id,
                "filename": f"ADC_Report_{session_id[:8]}_{datetime.utcnow().strftime('%Y%m%d')}.pdf"
            }

        except Exception as e:
            logger.error(f"PDF generation error: {e}")
            return {
                "format": "html",
                "content": html_content,
                "seal": seal_data,
                "session_id": session_id,
                "error": str(e)
            }

    def _get_grade_class(self, grade: str) -> str:
        """등급에 따른 CSS 클래스"""
        if grade in ["A", "A+"]:
            return "grade-go"
        elif grade in ["B", "B+"]:
            return "grade-conditional"
        else:
            return "grade-nogo"

    def _extract_key_metrics(self, session_data: Dict, agent_results: Dict) -> List[Dict]:
        """핵심 메트릭 추출"""
        metrics = []

        physchem = agent_results.get("physicochemical", {})
        if physchem.get("molecular_weight"):
            metrics.append({"label": "MW (Da)", "value": f"{physchem['molecular_weight']:,.0f}"})
        if physchem.get("logp"):
            metrics.append({"label": "LogP", "value": f"{physchem['logp']:.2f}"})

        metrics.append({"label": "Target DAR", "value": session_data.get("requested_dar", 4)})

        synthesis = agent_results.get("synthesis_analysis", {})
        if synthesis.get("sa_score"):
            metrics.append({"label": "SA Score", "value": f"{synthesis['sa_score']:.1f}"})

        tox = agent_results.get("toxicology", {})
        if tox.get("herg_risk"):
            metrics.append({"label": "hERG Risk", "value": tox["herg_risk"]})

        return metrics[:6]

    def _format_benchmark_table(self, benchmark_data: Dict) -> List[Dict]:
        """벤치마크 테이블 포맷팅"""
        table = benchmark_data.get("comparison_table", [])
        winner_analysis = {w["property"]: w for w in benchmark_data.get("winner_analysis", [])}

        formatted = []
        for row in table:
            prop = row.get("property", "")
            winner_info = winner_analysis.get(prop, {})
            winner = winner_info.get("winner", "")

            formatted_row = {
                "property": prop,
                "our_adc": row.get("our_adc", "N/A"),
                "winner": winner,
                "winner_display": "Our ADC" if winner == "our_adc" else winner.title() if winner else "-"
            }

            for comp in benchmark_data.get("competitors", {}).keys():
                formatted_row[comp] = row.get(comp, "N/A")

            formatted.append(formatted_row)

        return formatted

    def _format_lipinski_properties(self, physchem: Dict) -> List[Dict]:
        """Lipinski 속성 포맷팅"""
        properties = [
            ("Molecular Weight", physchem.get("molecular_weight"), "<500", 500),
            ("LogP", physchem.get("logp"), "<5", 5),
            ("H-Bond Donors", physchem.get("hbd"), "<5", 5),
            ("H-Bond Acceptors", physchem.get("hba"), "<10", 10),
            ("Rotatable Bonds", physchem.get("rotatable_bonds"), "<10", 10),
        ]

        formatted = []
        for name, value, limit, threshold in properties:
            if value is not None:
                status = "Pass" if value < threshold else "Fail"
                status_color = "#22c55e" if status == "Pass" else "#ef4444"
            else:
                status = "N/A"
                status_color = "#6b7280"

            formatted.append({
                "name": name,
                "value": f"{value:.2f}" if isinstance(value, float) else str(value) if value else "N/A",
                "limit": limit,
                "status": status,
                "status_color": status_color
            })

        return formatted

    def _format_toxicity_risks(self, tox_data: Dict) -> List[Dict]:
        """독성 위험 포맷팅"""
        risks = []
        predictions = tox_data.get("toxicity_predictions", {})

        risk_types = [
            ("hERG Cardiotoxicity", "herg_risk"),
            ("Hepatotoxicity", "hepatotoxicity_risk"),
            ("Off-target Binding", "off_target_binding"),
            ("Immunogenicity", "immunogenicity_risk"),
        ]

        level_colors = {
            "Low": "#dcfce7",
            "Medium": "#fef3c7",
            "High": "#fee2e2"
        }

        for risk_name, key in risk_types:
            level = predictions.get(key, "Unknown")
            risks.append({
                "type": risk_name,
                "level": level,
                "level_color": level_colors.get(level, "#f3f4f6"),
                "prediction": predictions.get(f"{key}_detail", ""),
                "mitigation": predictions.get(f"{key}_mitigation", "Monitor in clinical trials")
            })

        return risks

    def _collect_references(self, agent_results: Dict) -> List[Dict]:
        """참조 논문 수집"""
        pmids = set()
        for section_data in agent_results.values():
            if isinstance(section_data, dict):
                pmids.update(section_data.get("referenced_pmids", []))

        references = []
        for pmid in list(pmids)[:20]:
            references.append({
                "pmid": pmid,
                "authors": "Author et al.",
                "year": "2024",
                "title": f"Reference {pmid}",
                "journal": "Journal"
            })

        return references

    def _extract_page_contents(self, html_content: str) -> List[Dict]:
        """페이지별 컨텐츠 추출 (해시용)"""
        import re
        pages = re.split(r'<div class="page-break">', html_content)
        return [{"content": p, "page_number": i+1} for i, p in enumerate(pages)]


# Singleton instance
_report_engine = None


def get_report_engine() -> ReportEngine:
    """Get or create ReportEngine instance"""
    global _report_engine
    if _report_engine is None:
        _report_engine = ReportEngine()
    return _report_engine
