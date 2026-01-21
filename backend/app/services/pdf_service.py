"""
PDF Report Service - WeasyPrint 기반 PDF 생성
ADC 분석 결과를 전문 리포트로 변환
"""
import os
from typing import Dict, Any, Optional
from datetime import datetime
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML, CSS
import tempfile

from app.core.config import settings


# 템플릿 디렉토리
TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "templates")


def get_grade_color(grade: str) -> str:
    """등급별 색상 반환"""
    colors = {
        "A": "#22c55e",   # Green
        "B+": "#84cc16",  # Lime
        "B": "#eab308",   # Yellow
        "C": "#f97316",   # Orange
        "D": "#ef4444",   # Red
    }
    return colors.get(grade, "#6b7280")


def get_recommendation_badge(recommendation: str) -> Dict[str, str]:
    """권고사항 배지 스타일"""
    badges = {
        "Go": {"bg": "#dcfce7", "text": "#166534", "border": "#22c55e"},
        "Conditional Go": {"bg": "#fef9c3", "text": "#854d0e", "border": "#eab308"},
        "Proceed with Caution": {"bg": "#ffedd5", "text": "#9a3412", "border": "#f97316"},
        "No Go": {"bg": "#fee2e2", "text": "#991b1b", "border": "#ef4444"},
    }
    return badges.get(recommendation, {"bg": "#f3f4f6", "text": "#374151", "border": "#9ca3af"})


async def generate_pdf_report(
    job_id: str,
    result_data: Dict[str, Any],
    user_info: Optional[Dict[str, str]] = None
) -> str:
    """
    ADC 분석 결과를 PDF 리포트로 생성
    
    Args:
        job_id: 작업 ID
        result_data: 분석 결과 데이터 (ADCState에서 추출)
        user_info: 사용자 정보 (이름, 소속 등)
    
    Returns:
        생성된 PDF 파일 경로
    """
    # Jinja2 환경 설정
    env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = env.get_template("report_template.jinja2")
    
    # 템플릿 데이터 준비
    grade = result_data.get("final_grade", "B")
    recommendation = result_data.get("recommendation", "Conditional Go")
    
    template_data = {
        # 메타 정보
        "job_id": job_id,
        "generated_at": datetime.now().strftime("%Y년 %m월 %d일 %H:%M"),
        "user_name": user_info.get("name", "Anonymous") if user_info else "Anonymous",
        "user_org": user_info.get("organization", "") if user_info else "",
        
        # ADC 구성
        "input": result_data.get("input", {}),
        
        # 결과 요약
        "grade": grade,
        "grade_color": get_grade_color(grade),
        "recommendation": recommendation,
        "recommendation_badge": get_recommendation_badge(recommendation),
        "executive_summary": result_data.get("executive_summary", ""),
        
        # 점수
        "scores": result_data.get("scores", {}),
        
        # 상세 분석
        "structure_analysis": result_data.get("structure_analysis"),
        "toxicity_risks": result_data.get("toxicity_risks", []),
        "patent_landscape": result_data.get("patent_landscape", []),
        "competitors": result_data.get("competitors", []),
        "clinical_protocol": result_data.get("clinical_protocol"),
        
        # 스타일링
        "primary_color": "#6366f1",  # Indigo
        "accent_color": "#ec4899",   # Pink
    }
    
    # HTML 렌더링
    html_content = template.render(**template_data)
    
    # CSS 스타일
    css = CSS(string=get_report_css())
    
    # PDF 생성
    with tempfile.NamedTemporaryFile(suffix=".pdf", delete=False) as tmp:
        HTML(string=html_content).write_pdf(tmp.name, stylesheets=[css])
        return tmp.name


def get_report_css() -> str:
    """PDF 리포트용 CSS 스타일"""
    return """
    @page {
        size: A4;
        margin: 2cm;
        @top-right {
            content: "ADC-GenAI Report";
            font-size: 10px;
            color: #9ca3af;
        }
        @bottom-center {
            content: counter(page) " / " counter(pages);
            font-size: 10px;
            color: #9ca3af;
        }
    }
    
    body {
        font-family: 'Noto Sans KR', 'Inter', sans-serif;
        font-size: 11pt;
        line-height: 1.6;
        color: #1f2937;
    }
    
    h1 {
        font-size: 24pt;
        font-weight: 700;
        color: #111827;
        margin-bottom: 0.5em;
    }
    
    h2 {
        font-size: 16pt;
        font-weight: 600;
        color: #374151;
        border-bottom: 2px solid #e5e7eb;
        padding-bottom: 0.3em;
        margin-top: 1.5em;
    }
    
    h3 {
        font-size: 13pt;
        font-weight: 600;
        color: #4b5563;
        margin-top: 1em;
    }
    
    .header {
        text-align: center;
        margin-bottom: 2em;
    }
    
    .logo {
        font-size: 14pt;
        font-weight: 700;
        color: #6366f1;
    }
    
    .grade-box {
        display: inline-block;
        width: 80px;
        height: 80px;
        border-radius: 50%;
        text-align: center;
        line-height: 80px;
        font-size: 36pt;
        font-weight: 700;
        color: white;
    }
    
    .badge {
        display: inline-block;
        padding: 4px 12px;
        border-radius: 20px;
        font-size: 10pt;
        font-weight: 500;
    }
    
    .summary-box {
        background: #f9fafb;
        border: 1px solid #e5e7eb;
        border-radius: 8px;
        padding: 1em;
        margin: 1em 0;
    }
    
    .score-bar {
        height: 8px;
        background: #e5e7eb;
        border-radius: 4px;
        overflow: hidden;
    }
    
    .score-fill {
        height: 100%;
        border-radius: 4px;
    }
    
    table {
        width: 100%;
        border-collapse: collapse;
        margin: 1em 0;
    }
    
    th, td {
        padding: 8px 12px;
        text-align: left;
        border-bottom: 1px solid #e5e7eb;
    }
    
    th {
        background: #f3f4f6;
        font-weight: 600;
        font-size: 10pt;
    }
    
    .risk-high { background: #fee2e2; color: #991b1b; }
    .risk-medium { background: #fef3c7; color: #92400e; }
    .risk-low { background: #dcfce7; color: #166534; }
    
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
        font-size: 9pt;
        color: #92400e;
    }
    """
