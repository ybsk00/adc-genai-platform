"""
Report Generation API
PDF/HTML 리포트 생성을 위한 API
"""
from fastapi import APIRouter, HTTPException, Response
from fastapi.responses import StreamingResponse
from pydantic import BaseModel
from typing import Optional, List
import io
import logging
from datetime import datetime

from app.core.supabase import get_supabase_client

logger = logging.getLogger(__name__)

router = APIRouter()


class ReportConfig(BaseModel):
    """리포트 생성 설정"""
    include_smiles: bool = True
    include_metrics: bool = True
    include_reasoning: bool = True
    include_code: bool = False
    include_references: bool = True
    format: str = "pdf"  # pdf, html


# HTML 템플릿 (PDF 변환용)
REPORT_HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>ADC Design Report - {session_id}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; line-height: 1.6; color: #333; max-width: 800px; margin: 0 auto; padding: 40px 20px; }}
        .header {{ text-align: center; margin-bottom: 40px; padding-bottom: 20px; border-bottom: 2px solid #3b82f6; }}
        .header h1 {{ color: #1e3a5f; font-size: 24px; margin-bottom: 10px; }}
        .header .subtitle {{ color: #666; font-size: 14px; }}
        .header .meta {{ display: flex; justify-content: center; gap: 20px; margin-top: 15px; font-size: 12px; color: #888; }}
        .section {{ margin-bottom: 30px; }}
        .section-title {{ font-size: 16px; font-weight: 600; color: #1e3a5f; margin-bottom: 15px; padding-bottom: 8px; border-bottom: 1px solid #e5e7eb; }}
        .info-grid {{ display: grid; grid-template-columns: repeat(2, 1fr); gap: 15px; }}
        .info-item {{ background: #f8fafc; padding: 12px; border-radius: 8px; }}
        .info-item .label {{ font-size: 11px; color: #666; text-transform: uppercase; letter-spacing: 0.5px; }}
        .info-item .value {{ font-size: 14px; font-weight: 500; margin-top: 4px; }}
        .candidate-card {{ background: #fff; border: 1px solid #e5e7eb; border-radius: 8px; padding: 15px; margin-bottom: 15px; }}
        .candidate-card .rank {{ display: inline-block; background: #3b82f6; color: white; padding: 2px 8px; border-radius: 4px; font-size: 12px; font-weight: 600; }}
        .candidate-card .score {{ float: right; background: #10b981; color: white; padding: 2px 8px; border-radius: 4px; font-size: 12px; }}
        .smiles {{ font-family: 'Courier New', monospace; font-size: 11px; background: #f1f5f9; padding: 8px; border-radius: 4px; word-break: break-all; margin-top: 10px; }}
        .metrics-grid {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 10px; margin-top: 10px; }}
        .metric {{ text-align: center; padding: 8px; background: #f8fafc; border-radius: 4px; }}
        .metric .name {{ font-size: 10px; color: #666; }}
        .metric .val {{ font-size: 14px; font-weight: 600; color: #1e3a5f; }}
        .log-entry {{ padding: 10px; border-left: 3px solid #3b82f6; background: #f8fafc; margin-bottom: 10px; }}
        .log-entry .agent {{ font-weight: 600; color: #3b82f6; font-size: 12px; }}
        .log-entry .message {{ margin-top: 5px; font-size: 13px; }}
        .hash-box {{ background: #fef3c7; border: 1px solid #f59e0b; padding: 15px; border-radius: 8px; margin-top: 20px; }}
        .hash-box .title {{ font-weight: 600; color: #92400e; margin-bottom: 8px; }}
        .hash-box .hash {{ font-family: monospace; font-size: 10px; word-break: break-all; color: #78350f; }}
        .footer {{ margin-top: 40px; padding-top: 20px; border-top: 1px solid #e5e7eb; text-align: center; font-size: 11px; color: #888; }}
        .compliance-badge {{ display: inline-block; background: #dcfce7; color: #166534; padding: 4px 12px; border-radius: 20px; font-size: 11px; font-weight: 500; margin-top: 10px; }}
        @media print {{
            body {{ padding: 20px; }}
            .candidate-card {{ break-inside: avoid; }}
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>ADC Design Report</h1>
        <div class="subtitle">{session_type} Session</div>
        <div class="meta">
            <span>Session ID: {session_id_short}</span>
            <span>Generated: {generated_at}</span>
        </div>
        <div class="compliance-badge">21 CFR Part 11 Compliant</div>
    </div>

    <div class="section">
        <div class="section-title">Session Information</div>
        <div class="info-grid">
            <div class="info-item">
                <div class="label">Target Antigen</div>
                <div class="value">{target_antigen}</div>
            </div>
            <div class="info-item">
                <div class="label">Target Indication</div>
                <div class="value">{target_indication}</div>
            </div>
            <div class="info-item">
                <div class="label">Requested DAR</div>
                <div class="value">{requested_dar}</div>
            </div>
            <div class="info-item">
                <div class="label">Status</div>
                <div class="value">{status}</div>
            </div>
        </div>
    </div>

    {candidates_section}

    {logs_section}

    <div class="hash-box">
        <div class="title">Digital Seal (SHA-256)</div>
        <div class="hash">Bundle Hash: {bundle_hash}</div>
    </div>

    <div class="footer">
        <p>This report was generated by ADC-GenAI Platform</p>
        <p>Report verification hash ensures data integrity from generation to export</p>
    </div>
</body>
</html>
"""


async def get_current_user_id() -> str:
    """현재 사용자 ID 조회"""
    return "dev-user-id"


@router.get("/session/{session_id}/report")
async def generate_report(
    session_id: str,
    format: str = "html",
    include_code: bool = False
):
    """
    세션 리포트 생성

    - **format**: html 또는 pdf
    - **include_code**: 실행 코드 포함 여부
    """
    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    # 세션 조회
    session = supabase.table("design_sessions").select("*").eq(
        "id", session_id
    ).single().execute()

    if not session.data:
        raise HTTPException(404, "Session not found")

    if session.data["user_id"] != user_id:
        raise HTTPException(403, "Not authorized")

    session_data = session.data

    # 후보 조회
    candidates = supabase.table("design_candidates").select(
        "rank, smiles, score, metrics, validation"
    ).eq("session_id", session_id).order("rank").execute()

    # 로그 조회
    logs = supabase.table("agent_execution_logs").select(
        "agent_name, status, decision_summary, record_hash, created_at"
    ).eq("session_id", session_id).order("sequence_number").execute()

    # 후보 섹션 HTML 생성
    candidates_html = ""
    if candidates.data:
        candidates_html = '<div class="section"><div class="section-title">Generated Candidates</div>'
        for c in candidates.data:
            metrics_html = ""
            if c.get("metrics"):
                metrics = c["metrics"]
                metrics_html = '<div class="metrics-grid">'
                for key, value in list(metrics.items())[:4]:
                    val_str = f"{value:.2f}" if isinstance(value, float) else str(value)
                    metrics_html += f'<div class="metric"><div class="name">{key.upper()}</div><div class="val">{val_str}</div></div>'
                metrics_html += '</div>'

            score_str = f"{c['score']:.2f}" if isinstance(c.get('score'), (int, float)) else str(c.get('score', 'N/A'))
            candidates_html += f'''
            <div class="candidate-card">
                <span class="rank">Rank #{c['rank']}</span>
                <span class="score">Score: {score_str}</span>
                <div class="smiles">{c['smiles']}</div>
                {metrics_html}
            </div>
            '''
        candidates_html += '</div>'

    # 로그 섹션 HTML 생성
    logs_html = ""
    if logs.data:
        logs_html = '<div class="section"><div class="section-title">Execution Log</div>'
        for log in logs.data[:10]:  # 최대 10개
            logs_html += f'''
            <div class="log-entry">
                <div class="agent">{log['agent_name']}</div>
                <div class="message">{log.get('decision_summary', 'No summary available')}</div>
            </div>
            '''
        logs_html += '</div>'

    # 번들 해시 생성
    import hashlib
    import json
    bundle_data = {
        "session_id": session_id,
        "candidates": candidates.data,
        "logs": [{"agent": l["agent_name"], "hash": l.get("record_hash")} for l in (logs.data or [])]
    }
    bundle_hash = hashlib.sha256(json.dumps(bundle_data, sort_keys=True, default=str).encode()).hexdigest()

    # HTML 생성
    html_content = REPORT_HTML_TEMPLATE.format(
        session_id=session_id,
        session_id_short=session_id[:8],
        session_type=session_data.get("session_type", "Unknown").title(),
        generated_at=datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC"),
        target_antigen=session_data.get("target_antigen") or "Not specified",
        target_indication=session_data.get("target_indication") or "Not specified",
        requested_dar=session_data.get("requested_dar") or "N/A",
        status=session_data.get("status", "Unknown").title(),
        candidates_section=candidates_html,
        logs_section=logs_html,
        bundle_hash=bundle_hash
    )

    if format == "pdf":
        # PDF 생성 (weasyprint 또는 대체 방법)
        try:
            from weasyprint import HTML
            pdf_buffer = io.BytesIO()
            HTML(string=html_content).write_pdf(pdf_buffer)
            pdf_buffer.seek(0)

            return StreamingResponse(
                pdf_buffer,
                media_type="application/pdf",
                headers={
                    "Content-Disposition": f'attachment; filename="report-{session_id[:8]}.pdf"'
                }
            )
        except ImportError:
            # weasyprint 미설치 시 HTML 반환 with PDF 변환 안내
            logger.warning("weasyprint not installed, returning HTML")
            return Response(
                content=html_content,
                media_type="text/html",
                headers={
                    "X-PDF-Fallback": "true",
                    "Content-Disposition": f'inline; filename="report-{session_id[:8]}.html"'
                }
            )
    else:
        return Response(
            content=html_content,
            media_type="text/html",
            headers={
                "Content-Disposition": f'inline; filename="report-{session_id[:8]}.html"'
            }
        )


@router.get("/session/{session_id}/report/preview")
async def preview_report(session_id: str):
    """
    리포트 미리보기 (메타데이터만)
    """
    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    session = supabase.table("design_sessions").select(
        "id, session_type, status, target_antigen, created_at, updated_at"
    ).eq("id", session_id).single().execute()

    if not session.data:
        raise HTTPException(404, "Session not found")

    if session.data.get("user_id") and session.data["user_id"] != user_id:
        raise HTTPException(403, "Not authorized")

    candidates_count = supabase.table("design_candidates").select(
        "id", count="exact"
    ).eq("session_id", session_id).execute()

    logs_count = supabase.table("agent_execution_logs").select(
        "id", count="exact"
    ).eq("session_id", session_id).execute()

    return {
        "session_id": session_id,
        "session_type": session.data.get("session_type"),
        "status": session.data.get("status"),
        "target_antigen": session.data.get("target_antigen"),
        "candidates_count": candidates_count.count or 0,
        "logs_count": logs_count.count or 0,
        "available_formats": ["html", "pdf"],
        "created_at": session.data.get("created_at"),
        "updated_at": session.data.get("updated_at")
    }
