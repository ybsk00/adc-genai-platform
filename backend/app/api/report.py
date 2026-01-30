"""
Report Generation API
10페이지 정밀 보고서 생성 + Digital Seal + QR 검증
"""
from fastapi import APIRouter, HTTPException, Response, BackgroundTasks
from fastapi.responses import StreamingResponse, JSONResponse
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
import io
import logging
from datetime import datetime

from app.core.supabase import get_supabase_client

logger = logging.getLogger(__name__)


class ReportGenerateRequest(BaseModel):
    """보고서 생성 요청"""
    format: str = "pdf"
    include_appendix: bool = True
    include_benchmarks: bool = True
    include_visualizations: bool = True
    language: str = "en"

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


# ============================================
# Enhanced Report Engine Endpoints (v2.0)
# ============================================

@router.post("/generate/{session_id}")
async def generate_enhanced_report(
    session_id: str,
    request: ReportGenerateRequest,
    background_tasks: BackgroundTasks
):
    """
    10페이지 정밀 보고서 생성 (Enhanced v2.0)

    Features:
    - AI Reasoning Logic 포함
    - Confidence Score (0-100%)
    - Digital Seal + QR 코드
    - 벤치마크 비교표
    - 시각화 자동 생성
    """
    user_id = await get_current_user_id()
    supabase = get_supabase_client()

    # 세션 조회
    session = supabase.table("design_sessions").select("*").eq(
        "id", session_id
    ).single().execute()

    if not session.data:
        raise HTTPException(404, "Session not found")

    session_data = session.data

    # 에이전트 실행 결과 조회
    logs = supabase.table("agent_execution_logs").select("*").eq(
        "session_id", session_id
    ).order("sequence_number").execute()

    # 후보 조회
    candidates = supabase.table("design_candidates").select("*").eq(
        "session_id", session_id
    ).order("rank").execute()

    # 에이전트 결과 집계
    agent_results = _aggregate_agent_results(logs.data or [], candidates.data or [])

    try:
        from app.services.report_engine import get_report_engine

        engine = get_report_engine()
        result = await engine.generate_report(
            session_id=session_id,
            session_data=session_data,
            agent_results=agent_results,
            format=request.format,
            include_visualizations=request.include_visualizations,
            include_benchmarks=request.include_benchmarks
        )

        # 감사 로그 저장
        background_tasks.add_task(
            _save_report_audit_log,
            session_id, user_id, result.get("seal", {})
        )

        if result["format"] == "pdf":
            return StreamingResponse(
                io.BytesIO(result["content"]),
                media_type="application/pdf",
                headers={
                    "Content-Disposition": f'attachment; filename="{result["filename"]}"',
                    "X-Report-Hash": result["seal"]["chain_hash"][:12],
                    "X-Verification-URL": result["seal"]["verification_url"]
                }
            )
        else:
            return Response(
                content=result["content"],
                media_type="text/html",
                headers={
                    "Content-Disposition": f'inline; filename="report-{session_id[:8]}.html"',
                    "X-Report-Hash": result["seal"]["chain_hash"][:12] if result.get("seal") else ""
                }
            )

    except Exception as e:
        logger.error(f"Report generation error: {e}")
        raise HTTPException(500, f"Report generation failed: {str(e)}")


@router.get("/verify/{hash}")
async def verify_report(hash: str):
    """
    보고서 무결성 검증 (QR 코드 스캔용)

    Args:
        hash: 보고서 해시 (12자리 short hash)

    Returns:
        검증 결과 및 보고서 메타데이터
    """
    supabase = get_supabase_client()

    # 보고서 해시로 검색
    report_log = supabase.table("report_verification_logs").select("*").eq(
        "short_hash", hash
    ).single().execute()

    if not report_log.data:
        # short_hash가 chain_hash의 앞 12자리인 경우도 검색
        report_log = supabase.table("report_verification_logs").select("*").ilike(
            "chain_hash", f"{hash}%"
        ).single().execute()

    if not report_log.data:
        return JSONResponse(
            status_code=404,
            content={
                "valid": False,
                "error": "Report not found",
                "hash": hash,
                "message": "This report hash could not be verified. The report may have been modified or does not exist in our records."
            }
        )

    # 세션 정보 조회
    session = supabase.table("design_sessions").select(
        "id, session_type, target_antigen, status, created_at"
    ).eq("id", report_log.data["session_id"]).single().execute()

    return {
        "valid": True,
        "report_id": report_log.data.get("report_id"),
        "session_id": report_log.data["session_id"],
        "generated_at": report_log.data["generated_at"],
        "page_count": report_log.data.get("page_count", 10),
        "chain_hash": report_log.data["chain_hash"],
        "algorithm": "SHA-256 Chain",
        "verification_timestamp": datetime.utcnow().isoformat(),
        "session_info": {
            "type": session.data.get("session_type") if session.data else None,
            "target": session.data.get("target_antigen") if session.data else None,
            "status": session.data.get("status") if session.data else None
        },
        "compliance": {
            "standard": "21 CFR Part 11",
            "features": ["Digital Seal", "Audit Trail", "QR Verification"]
        }
    }


@router.get("/verify/{hash}/details")
async def get_verification_details(hash: str):
    """
    보고서 검증 상세 정보

    각 페이지별 해시 및 감사 추적 정보 반환
    """
    supabase = get_supabase_client()

    report_log = supabase.table("report_verification_logs").select("*").ilike(
        "chain_hash", f"{hash}%"
    ).single().execute()

    if not report_log.data:
        raise HTTPException(404, "Report not found")

    return {
        "chain_hash": report_log.data["chain_hash"],
        "page_hashes": report_log.data.get("page_hashes", []),
        "generated_at": report_log.data["generated_at"],
        "generated_by": report_log.data.get("generated_by"),
        "audit_entries": report_log.data.get("audit_trail", [])
    }


def _aggregate_agent_results(logs: List[Dict], candidates: List[Dict]) -> Dict[str, Any]:
    """에이전트 결과 집계"""
    results = {
        "target_analysis": {},
        "synthesis_analysis": {},
        "benchmark_comparison": {},
        "physicochemical": {},
        "toxicology": {},
        "patent_analysis": {},
        "audit_trail": [],
        "final_grade": "B",
        "recommendation": "Conditional Go",
        "executive_summary": "",
        "executive_reasoning": ""
    }

    for log in logs:
        agent_name = log.get("agent_name", "").lower()
        output = log.get("output_data", {}) or {}

        # 에이전트별 결과 매핑
        if "librarian" in agent_name:
            results["target_analysis"] = output
        elif "alchemist" in agent_name:
            results["synthesis_analysis"] = output
        elif "coder" in agent_name:
            results["physicochemical"] = output
        elif "healer" in agent_name or "tox" in agent_name:
            results["toxicology"] = output
        elif "competitor" in agent_name or "benchmark" in agent_name:
            results["benchmark_comparison"] = output
        elif "patent" in agent_name:
            results["patent_analysis"] = output
        elif "auditor" in agent_name:
            results["audit_trail"] = output.get("audit_entries", [])

        # 감사 추적 항목 추가
        results["audit_trail"].append({
            "timestamp": log.get("created_at"),
            "action": log.get("status"),
            "agent": log.get("agent_name"),
            "hash": log.get("record_hash", "")[:12] if log.get("record_hash") else ""
        })

    # 최종 판정 (오케스트레이터 결과 또는 기본값)
    orchestrator_log = next(
        (l for l in logs if "orchestrator" in l.get("agent_name", "").lower()),
        None
    )
    if orchestrator_log:
        output = orchestrator_log.get("output_data", {}) or {}
        results["final_grade"] = output.get("final_grade", "B")
        results["recommendation"] = output.get("recommendation", "Conditional Go")
        results["executive_summary"] = output.get("summary", "")
        results["executive_reasoning"] = output.get("reasoning_logic", "")

    # 후보 데이터에서 물리화학적 속성 추출
    if candidates:
        top_candidate = candidates[0]
        metrics = top_candidate.get("metrics", {}) or {}
        results["physicochemical"].update({
            "molecular_weight": metrics.get("mw") or metrics.get("molecular_weight"),
            "logp": metrics.get("logp") or metrics.get("LogP"),
            "hbd": metrics.get("hbd"),
            "hba": metrics.get("hba"),
            "tpsa": metrics.get("tpsa"),
            "rotatable_bonds": metrics.get("rotatable_bonds")
        })

    return results


async def _save_report_audit_log(
    session_id: str,
    user_id: str,
    seal_data: Dict
):
    """보고서 감사 로그 저장"""
    try:
        supabase = get_supabase_client()

        supabase.table("report_verification_logs").insert({
            "session_id": session_id,
            "report_id": f"RPT-{session_id[:8]}-{datetime.utcnow().strftime('%Y%m%d%H%M%S')}",
            "chain_hash": seal_data.get("chain_hash", ""),
            "short_hash": seal_data.get("short_hash", ""),
            "page_hashes": seal_data.get("page_hashes", []),
            "page_count": seal_data.get("page_count", 10),
            "generated_at": seal_data.get("generated_at"),
            "generated_by": user_id,
            "algorithm": seal_data.get("algorithm", "SHA-256 Chain"),
            "verification_url": seal_data.get("verification_url", "")
        }).execute()

        logger.info(f"Report audit log saved for session {session_id}")

    except Exception as e:
        logger.error(f"Failed to save report audit log: {e}")


# ============================================
# Phase 2: Digital Lineage API Endpoints
# ============================================

@router.get("/session/{session_id}/lineage")
async def get_session_lineage(session_id: str):
    """
    세션의 Digital Lineage (데이터 계보) 조회

    논문 작성 및 특허 출원 시 연산 환경 재현성 검증용
    """
    supabase = get_supabase_client()

    # 세션 확인
    session = supabase.table("design_sessions").select("id").eq(
        "id", session_id
    ).single().execute()

    if not session.data:
        raise HTTPException(404, "Session not found")

    # Calculation Lineage 조회
    lineage = supabase.table("agent_execution_logs").select(
        "calculation_id, agent_name, agent_version, status, "
        "execution_env, library_versions, nvidia_nim_info, "
        "simulation_params, input_hash, output_hash, "
        "llm_model, llm_temperature, execution_time_ms, "
        "created_at"
    ).eq("session_id", session_id).eq("status", "completed").order("created_at").execute()

    # Physical Validations 조회
    validations = supabase.table("physical_validations").select(
        "calculation_id, check_name, check_category, passed, severity, details"
    ).eq("session_id", session_id).execute()

    # 결과 포맷팅
    lineage_data = []
    for log in lineage.data or []:
        calc_id = log.get("calculation_id")
        calc_validations = [v for v in (validations.data or []) if v.get("calculation_id") == calc_id]

        lineage_data.append({
            "calculation_id": calc_id,
            "timestamp": log.get("created_at"),
            "agent_info": {
                "name": log.get("agent_name"),
                "version": log.get("agent_version"),
                "llm_model": log.get("llm_model"),
                "llm_temperature": log.get("llm_temperature")
            },
            "execution_environment": log.get("execution_env", {}),
            "library_versions": log.get("library_versions", {}),
            "nvidia_nim": log.get("nvidia_nim_info"),
            "simulation_params": log.get("simulation_params"),
            "hashes": {
                "input": log.get("input_hash"),
                "output": log.get("output_hash")
            },
            "execution_time_ms": log.get("execution_time_ms"),
            "validations": {
                "total": len(calc_validations),
                "passed": sum(1 for v in calc_validations if v.get("passed")),
                "failed": sum(1 for v in calc_validations if not v.get("passed")),
                "details": calc_validations[:5]  # 최대 5개만
            }
        })

    return {
        "session_id": session_id,
        "total_calculations": len(lineage_data),
        "lineage": lineage_data,
        "compliance": {
            "standard": "21 CFR Part 11",
            "features": ["Digital Lineage", "Input/Output Hash", "Library Versioning"]
        }
    }


@router.get("/lineage/{calculation_id}")
async def get_calculation_metadata(calculation_id: str):
    """
    특정 계산의 상세 메타데이터 조회

    논문/특허에서 특정 수치 인용 시 사용
    """
    supabase = get_supabase_client()

    # Calculation Metadata 조회
    log = supabase.table("agent_execution_logs").select("*").eq(
        "calculation_id", calculation_id
    ).single().execute()

    if not log.data:
        raise HTTPException(404, "Calculation not found")

    # Physical Validations 조회
    validations = supabase.table("physical_validations").select("*").eq(
        "calculation_id", calculation_id
    ).execute()

    return {
        "calculation_id": calculation_id,
        "timestamp": log.data.get("created_at"),
        "execution_environment": {
            "sandbox_type": log.data.get("execution_env", {}).get("sandbox_type", "subprocess"),
            "container_image": log.data.get("execution_env", {}).get("container_image", "N/A"),
            "python_version": log.data.get("execution_env", {}).get("python_version", "N/A"),
            "os_info": log.data.get("execution_env", {}).get("os_info", "N/A")
        },
        "library_versions": log.data.get("library_versions", {}),
        "nvidia_nim": log.data.get("nvidia_nim_info"),
        "simulation_params": log.data.get("simulation_params"),
        "agent_info": {
            "agent_name": log.data.get("agent_name"),
            "agent_version": log.data.get("agent_version"),
            "llm_model": log.data.get("llm_model"),
            "llm_temperature": log.data.get("llm_temperature")
        },
        "hashes": {
            "input": log.data.get("input_hash"),
            "output": log.data.get("output_hash")
        },
        "physical_validations": validations.data or [],
        "reproducibility_info": {
            "note": "To reproduce this calculation, use the same library versions and input data",
            "input_hash": log.data.get("input_hash"),
            "expected_output_hash": log.data.get("output_hash")
        }
    }


class ReproducibilityVerifyRequest(BaseModel):
    """재현성 검증 요청"""
    calculation_id: str
    new_output_hash: str


@router.post("/lineage/verify-reproducibility")
async def verify_reproducibility(request: ReproducibilityVerifyRequest):
    """
    계산 재현성 검증

    동일한 입력으로 다시 계산한 결과의 output_hash가 일치하는지 확인
    """
    supabase = get_supabase_client()

    # Original calculation 조회
    log = supabase.table("agent_execution_logs").select(
        "calculation_id, output_hash, execution_env, library_versions, created_at"
    ).eq("calculation_id", request.calculation_id).single().execute()

    if not log.data:
        raise HTTPException(404, "Original calculation not found")

    original_hash = log.data.get("output_hash")
    is_reproducible = original_hash == request.new_output_hash

    # 검증 로그 저장
    try:
        supabase.table("report_verification_logs").insert({
            "report_id": f"REPRO-{request.calculation_id[:8]}",
            "session_id": log.data.get("session_id"),
            "verification_type": "reproducibility",
            "verification_result": "valid" if is_reproducible else "invalid",
            "chain_hash_matched": is_reproducible,
            "details": {
                "original_hash": original_hash,
                "new_hash": request.new_output_hash,
                "execution_env": log.data.get("execution_env"),
                "library_versions": log.data.get("library_versions")
            },
            "verified_at": datetime.utcnow().isoformat()
        }).execute()
    except Exception as e:
        logger.warning(f"Failed to save reproducibility log: {e}")

    return {
        "calculation_id": request.calculation_id,
        "is_reproducible": is_reproducible,
        "original_hash": original_hash,
        "new_hash": request.new_output_hash,
        "original_timestamp": log.data.get("created_at"),
        "verification_timestamp": datetime.utcnow().isoformat(),
        "execution_environment": log.data.get("execution_env", {}),
        "library_versions": log.data.get("library_versions", {}),
        "message": "Calculation is reproducible!" if is_reproducible else "Hash mismatch - calculation may have been affected by different library versions or environment"
    }
