"""
Data Worker API - 외부 데이터 소스 동기화 (AstraForge 2.0)
DB 기반 상태 관리(Supabase) 및 PubChem 화학 정보 연동
"""
from fastapi import APIRouter, BackgroundTasks, HTTPException
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
from datetime import datetime, timedelta, timezone
from uuid import uuid4
import aiohttp
import asyncio
from Bio import Entrez
import logging
import json
import subprocess
import sys
import os

from app.core.config import settings
from app.core.supabase import supabase
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from app.services.openfda_service import openfda_service
from app.services.scheduler_engine import scheduler_engine
from app.services.chemical_resolver import chemical_resolver
from app.services.job_lock import job_lock

router = APIRouter()
logger = logging.getLogger(__name__)

# --- Job Manager (In-Memory Tracking) ---
class JobManager:
    def __init__(self):
        self.active_jobs = set()

    def add_job(self, job_id: str):
        self.active_jobs.add(job_id)

    def remove_job(self, job_id: str):
        if job_id in self.active_jobs:
            self.active_jobs.remove(job_id)

    def is_active(self, job_id: str) -> bool:
        return job_id in self.active_jobs

job_manager = JobManager()

# Entrez 설정
Entrez.email = settings.NCBI_EMAIL
if settings.NCBI_API_KEY:
    Entrez.api_key = settings.NCBI_API_KEY

class DataSourceSetting(BaseModel):
    source_id: str
    auto_sync: bool
    sync_interval_hours: int

class SyncJobResponse(BaseModel):
    job_id: str
    status: str
    message: str

class SyncJobStatus(BaseModel):
    id: str
    status: str
    source: str
    records_found: int
    records_drafted: int
    last_processed_page: Optional[int] = 0
    started_at: str
    completed_at: Optional[str] = None
    errors: List[str] = []

# --- DB 기반 상태 관리 유틸리티 ---

async def update_job_status(job_id: str, **kwargs):
    """DB의 sync_jobs 테이블 상태 업데이트"""
    try:
        # errors 필드가 있으면 JSON 직렬화 확인
        if "errors" in kwargs and isinstance(kwargs["errors"], list):
            kwargs["errors"] = kwargs["errors"] # Supabase client handles list to jsonb
        
        # message 필드는 sync_jobs 테이블에 없으므로 제거
        if "message" in kwargs:
            del kwargs["message"]

        # DB 업데이트 실행
        try:
            supabase.table("sync_jobs").update(kwargs).eq("id", job_id).execute()
        except Exception as e:
            # last_processed_page 컬럼이 없는 경우를 대비한 자동 폴백
            if "last_processed_page" in str(e) or "42703" in str(e):
                logger.warning("⚠️ 'last_processed_page' column missing. Retrying without it.")
                if "last_processed_page" in kwargs:
                    del kwargs["last_processed_page"]
                supabase.table("sync_jobs").update(kwargs).eq("id", job_id).execute()
            else:
                raise e
    except Exception as e:
        logger.error(f"Failed to update job status in DB: {e}")

async def get_job_from_db(job_id: str) -> Optional[dict]:
    """DB에서 작업 상태 조회"""
    try:
        res = supabase.table("sync_jobs").select("*").eq("id", job_id).execute()
        return res.data[0] if res.data else None
    except Exception as e:
        logger.error(f"Failed to get job from DB: {e}")
        return None

async def is_cancelled(job_id: str) -> bool:
    """중단 요청 여부 확인"""
    job = await get_job_from_db(job_id)
    return job.get("cancel_requested", False) if job else False

async def run_isolated_crawler(crawler_type: str, category: str, limit: int, job_id: str, start_page: int = 1):
    """실행 격리를 위해 별도 프로세스로 크롤러 실행"""
    log_file = f"crawler_debug_{job_id}.log"
    try:
        # 절대 경로 계산
        base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        script_path = os.path.join(base_dir, "run_crawler.py")
        
        if not os.path.exists(script_path):
            logger.error(f"❌ Crawler script not found at: {script_path}")
            await update_job_status(job_id, status="failed", errors=[f"Script not found: {script_path}"])
            return

        cmd = [
            sys.executable,
            script_path,
            "--crawler", crawler_type,
            "--category", category,
            "--limit", str(limit),
            "--job_id", job_id,
            "--start_page", str(start_page)
        ]
        
        logger.info(f"🚀 Launching isolated crawler (Absolute Path): {' '.join(cmd)}")
        
        # 로그 파일 생성
        with open(log_file, "w") as f:
            f.write(f"Started at: {datetime.now()}\nCommand: {' '.join(cmd)}\n\n")

        process = await asyncio.create_subprocess_exec(
            *cmd,
            stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            env={**os.environ, "GRPC_ENABLE_FORK_SUPPORT": "false", "PYTHONPATH": base_dir}
        )
        
        # 비동기로 로그 캡처 (서버 로그에 출력)
        stdout, stderr = await process.communicate()
        
        if stdout: logger.info(f"📝 Crawler Output: {stdout.decode()[:500]}")
        if stderr: logger.error(f"⚠️ Crawler Error Output: {stderr.decode()[:500]}")

        if process.returncode != 0:
            error_msg = stderr.decode().strip() if stderr else "Unknown error"
            logger.error(f"❌ Crawler process failed (code {process.returncode}): {error_msg}")
            await update_job_status(job_id, status="failed", errors=[f"Exit code {process.returncode}", error_msg])
        else:
            logger.info(f"✅ Crawler process finished successfully.")

    except Exception as e:
        logger.error(f"🔥 Failed to launch isolated crawler: {e}", exc_info=True)
        await update_job_status(job_id, status="failed", errors=[str(e)])
    finally:
        job_manager.remove_job(job_id)

# --- LLM & Chemical Enrichment ---

async def refine_drug_data_with_llm(title: str, text_content: str, raw_status: Optional[str] = None) -> dict:
    """LLM을 사용하여 ADC 정보 추출 및 PubChem 보강"""
    try:
        llm_refiner = ChatOpenAI(
            model=settings.FAST_LLM,
            temperature=0,
            api_key=settings.OPENAI_API_KEY
        )
        
        system_prompt = """You are an expert Clinical & Pharmaceutical Data Analyst. 
        Extract structured ADC (Antibody-Drug Conjugate) information.
        Output JSON format:
        {
            "drug_name": {"value": "DS-8201a", "evidence": "..."},
            "target": {"value": "HER2", "evidence": "..."},
            "payload": {"value": "DXd", "evidence": "..."},
            "linker": {"value": "Tetrapeptide-based", "evidence": "..."},
            "outcome_type": {"value": "Success", "evidence": "..."},
            "smiles_code": {"value": null, "evidence": null}
        }"""
        
        user_prompt = f"Title: {title}\nContent: {text_content}\nStatus: {raw_status}"
        prompt = ChatPromptTemplate.from_messages([("system", system_prompt), ("user", user_prompt)])
        chain = prompt | llm_refiner
        response = await chain.ainvoke({})
        content = response.content.strip()
        if "```json" in content: content = content.split("```json")[1].split("```")[0]
        data = json.loads(content)
        
        # SMILES 보강 (PubChem Resolver)
        drug_name = data.get("drug_name", {}).get("value")
        if drug_name and not data.get("smiles_code", {}).get("value"):
            smiles = chemical_resolver.fetch_verified_smiles(drug_name)
            if smiles:
                data["smiles_code"] = {"value": smiles, "evidence": "Verified via PubChem API"}
        
        return data
    except Exception as e:
        logger.error(f"LLM Refine Error: {e}")
        return {"drug_name": {"value": "Unknown"}}

# --- ClinicalTrials.gov Worker (DB Status & Paging) ---

async def process_clinical_trials_data(job_id: str, max_records: int = 1000):
    """ClinicalTrials.gov 대량 수집 워커 (BulkImporter 통합)"""
    from app.services.bulk_importer import BulkImporter
    
    importer = BulkImporter()
    # 기본적으로 daily 모드로 동작 (필요 시 mode 파라미터 확장 가능)
    await importer.run_import(job_id, max_studies=max_records, mode="daily")

async def process_single_study_simple(study: dict) -> bool:
    """단순화된 스터디 저장 (LLM 없이)"""
    try:
        protocol = study.get("protocolSection", {})
        id_module = protocol.get("identificationModule", {})
        nct_id = id_module.get("nctId")
        if not nct_id: return False
        
        existing = supabase.table("golden_set_library").select("id").eq("properties->>nct_id", nct_id).execute()
        if existing.data: return False

        title = id_module.get("officialTitle") or id_module.get("briefTitle", "No Title")
        status_module = protocol.get("statusModule", {})
        
        new_entry = {
            "name": title[:100],
            "category": "clinical_trial",
            "description": title,
            "properties": {"nct_id": nct_id, "phase": status_module.get("phase"), "status": status_module.get("overallStatus")},
            "status": "draft",
            "outcome_type": "Unknown",
            "enrichment_source": "clinical_trials_simple"
        }
        supabase.table("golden_set_library").insert(new_entry).execute()
        return True
    except Exception as e:
        logger.error(f"Study save error {study.get('protocolSection', {}).get('identificationModule', {}).get('nctId')}: {e}")
        return False

async def process_single_study(study: dict, job_id: str) -> bool:
    """LLM 포함 스터디 처리 (기존 로직 유지)"""
    try:
        protocol = study.get("protocolSection", {})
        id_module = protocol.get("identificationModule", {})
        nct_id = id_module.get("nctId")
        if not nct_id: return False
        
        existing = supabase.table("golden_set_library").select("id").eq("properties->>nct_id", nct_id).execute()
        if existing.data: return False

        refined = await refine_drug_data_with_llm(
            title=id_module.get("officialTitle", "No Title"),
            text_content=str(protocol.get("descriptionModule", {})),
            raw_status=protocol.get("statusModule", {}).get("overallStatus")
        )
        
        new_entry = {
            "name": refined.get("drug_name", {}).get("value", "Unknown"),
            "category": "clinical_trial",
            "description": id_module.get("officialTitle", "No Title"),
            "properties": {"nct_id": nct_id, "ai_analysis": refined},
            "smiles_code": refined.get("smiles_code", {}).get("value"),
            "status": "draft",
            "outcome_type": refined.get("outcome_type", {}).get("value", "Unknown"),
            "enrichment_source": "clinical_trials_v2_db_status"
        }
        supabase.table("golden_set_library").insert(new_entry).execute()
        return True
    except Exception:
        return False

# --- PubMed Worker (DB Status & Paging) ---

async def process_pubmed_data(job_id: str, max_records: int = 500, mode: str = "daily"):
    """PubMed 대량 수집 워커 (DB 기반)"""
    await update_job_status(job_id, status="running")
    drafted = 0
    errors = []
    
    # 모드에 따른 설정
    days_back = 2 if mode == "daily" else 3650 # Daily: 2일, Full: 10년 (기본값)
    
    # Smart Trawl 전략: 2015년 이후 데이터만 (Full Load 시)
    mindate = "2015/01/01" if mode == "full" else None
    
    if mode == "full":
        max_records = 2000 # Full load 시 더 많이 수집
    
    try:
        loop = asyncio.get_event_loop()
        
        # Smart Trawl Query Construction
        # 1. Target Keywords
        target_terms = '(Antibody-Drug Conjugate[Title/Abstract] OR ADC[Title/Abstract] OR Immunoconjugate[Title/Abstract])'
        # 2. Context Keywords (Safety Guard)
        context_terms = '(Cancer[Title/Abstract] OR Tumor[Title/Abstract] OR Oncology[Title/Abstract] OR Neoplasms[MeSH Terms])'
        # 3. Final Query
        search_term = f"{target_terms} AND {context_terms}"
        
        id_list = await loop.run_in_executor(None, lambda: fetch_pubmed_ids(search_term, max_records, days_back, mindate))
        await update_job_status(job_id, records_found=len(id_list))
        
        if not id_list:
            await update_job_status(job_id, status="completed", completed_at=datetime.utcnow().isoformat())
            return

        for i in range(0, len(id_list), 10):
            if await is_cancelled(job_id):
                await update_job_status(job_id, status="stopped")
                return

            batch_ids = id_list[i:i+10]
            papers = await loop.run_in_executor(None, fetch_pubmed_details, batch_ids)
            
            if 'PubmedArticle' in papers:
                tasks = [process_single_article(article, job_id) for article in papers['PubmedArticle']]
                results = await asyncio.gather(*tasks, return_exceptions=True)
                
                for r in results:
                    if isinstance(r, Exception): errors.append(str(r))
                    elif r: drafted += 1
                
                await update_job_status(job_id, records_drafted=drafted, errors=errors[:20])

        await update_job_status(job_id, status="completed", completed_at=datetime.utcnow().isoformat())
    except Exception as e:
        logger.error(f"PubMed Worker Error: {e}")
        await update_job_status(job_id, status="failed", errors=[str(e)])

async def process_single_article(article: dict, job_id: str) -> bool:
    """PubMed 논문을 knowledge_base에 저장 (RAG용 지식 베이스)"""
    try:
        medline = article['MedlineCitation']
        article_data = medline['Article']
        title = str(article_data.get('ArticleTitle', 'No Title'))
        
        # 중복 체크 (title 기준)
        existing = supabase.table("knowledge_base").select("id").eq("title", title).execute()
        if existing.data: 
            return False
        
        # 초록 추출
        abstract_texts = article_data.get('Abstract', {}).get('AbstractText', [])
        if isinstance(abstract_texts, list):
            content = " ".join([str(t) for t in abstract_texts])
        else:
            content = str(abstract_texts) if abstract_texts else ""
        
        # PMID 추출
        pmid = str(medline.get('PMID', ''))
        
        # 저널 정보
        journal = article_data.get('Journal', {}).get('Title', '')
        
        # knowledge_base 스키마에 맞게 저장
        new_kb = {
            "source_type": "PubMed",
            "title": title[:500],
            "summary": f"PMID: {pmid} | Journal: {journal}" if pmid else "",
            "content": content if content else "No abstract available.",
            "relevance_score": None,
            "source_tier": 1,  # PubMed = Tier 1 (학술 논문)
            "ai_reasoning": None,
            "rag_status": "pending"
        }
        
        supabase.table("knowledge_base").insert(new_kb).execute()
        logger.info(f"✅ Saved PubMed: {title[:50]}...")
        return True
    except Exception as e:
        logger.error(f"PubMed save error: {e}")
        return False

def fetch_pubmed_ids(term: str, max_results: int, days_back: int = 7, mindate: Optional[str] = None):
    # 날짜 필터링
    end_date = datetime.now()
    
    # mindate가 명시되면 그것을 사용, 아니면 days_back 사용
    if mindate:
        start_date_str = mindate
    else:
        start_date = end_date - timedelta(days=days_back)
        start_date_str = start_date.strftime("%Y/%m/%d")
    
    handle = Entrez.esearch(
        db="pubmed", 
        term=term, 
        retmax=max_results, 
        sort="pub_date",
        datetype="pdat",
        mindate=start_date_str,
        maxdate=end_date.strftime("%Y/%m/%d")
    )
    record = Entrez.read(handle)
    handle.close()
    return record.get("IdList", [])

def fetch_pubmed_details(id_list: List[str]):
    handle = Entrez.efetch(db="pubmed", id=",".join(id_list), retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

# --- API Endpoints ---

@router.post("/sync/clinical", response_model=SyncJobResponse)
async def sync_clinical_trials(background_tasks: BackgroundTasks):
    job_id = f"sync_clinical_{uuid4().hex[:8]}"
    data = {"id": job_id, "status": "queued", "source": "clinical_trials", "started_at": datetime.utcnow().isoformat()}
    supabase.table("sync_jobs").insert(data).execute()
    background_tasks.add_task(process_clinical_trials_data, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="ClinicalTrials.gov sync started.")

@router.post("/sync/pubmed", response_model=SyncJobResponse)
async def sync_pubmed(background_tasks: BackgroundTasks, mode: str = "daily"):
    job_id = f"sync_pubmed_{uuid4().hex[:8]}"
    data = {"id": job_id, "status": "queued", "source": "pubmed", "started_at": datetime.utcnow().isoformat()}
    supabase.table("sync_jobs").insert(data).execute()
    background_tasks.add_task(process_pubmed_data, job_id, 500, mode) # Pass mode
    return SyncJobResponse(job_id=job_id, status="queued", message=f"PubMed sync started (mode: {mode}).")

@router.post("/sync/openfda", response_model=SyncJobResponse)
async def sync_openfda(background_tasks: BackgroundTasks, mode: str = "daily", limit: int = 100):
    job_id = f"sync_openfda_{uuid4().hex[:8]}"
    data = {"id": job_id, "status": "queued", "source": "openfda", "started_at": datetime.utcnow().isoformat()}
    supabase.table("sync_jobs").insert(data).execute()
    background_tasks.add_task(openfda_service.sync_to_db, job_id, mode, limit)
    return SyncJobResponse(job_id=job_id, status="queued", message=f"OpenFDA sync started (mode: {mode}, limit: {limit}).")

@router.post("/crawler/creative-biolabs/run", response_model=SyncJobResponse)
async def run_creative_biolabs_crawler(background_tasks: BackgroundTasks, category: str = "ADC Cytotoxin", limit: int = 10):
    """Creative Biolabs Stealth Crawler 실행 (격리 프로세스 방식)"""
    job_id = f"crawl_cb_{uuid4().hex[:8]}"
    
    # DB에 작업 기록
    data = {
        "id": job_id, 
        "status": "running", 
        "source": "creative_biolabs", 
        "started_at": datetime.now(timezone.utc).isoformat()
    }
    supabase.table("sync_jobs").insert(data).execute()
    
    job_manager.add_job(job_id)
    background_tasks.add_task(run_isolated_crawler, "creative_biolabs", category, limit, job_id)
    
    return SyncJobResponse(job_id=job_id, status="queued", message="Creative Biolabs crawler started (Isolated Process).")

@router.post("/crawler/ambeed/run", response_model=SyncJobResponse)
async def run_ambeed_crawler(background_tasks: BackgroundTasks, category: str = "all", limit: int = 10, start_page: int = 1):
    """Ambeed Stealth Crawler 실행 (격리 프로세스 방식)"""
    job_id = f"crawl_ambeed_{uuid4().hex[:8]}"
    
    data = {
        "id": job_id, 
        "status": "running", 
        "source": "ambeed", 
        "started_at": datetime.now(timezone.utc).isoformat()
    }
    supabase.table("sync_jobs").insert(data).execute()
    
    job_manager.add_job(job_id)
    background_tasks.add_task(run_isolated_crawler, "ambeed", category, limit, job_id, start_page)
    
    return SyncJobResponse(job_id=job_id, status="queued", message=f"Ambeed crawler started (Isolated Process) from page {start_page}.")

@router.get("/sync/{job_id}")
async def get_sync_status(job_id: str):
    job = await get_job_from_db(job_id)
    if not job: raise HTTPException(status_code=404, detail="Job not found")
    
    # Zombie Check
    if job.get("status") == "running" and not job_manager.is_active(job_id):
        started_at_str = job.get("started_at")
        try:
            # Handle ISO format with 'Z' or offset
            if started_at_str.endswith('Z'):
                started_at = datetime.fromisoformat(started_at_str[:-1]).replace(tzinfo=timezone.utc)
            else:
                started_at = datetime.fromisoformat(started_at_str)
                if started_at.tzinfo is None:
                    started_at = started_at.replace(tzinfo=timezone.utc)

            if (datetime.now(timezone.utc) - started_at).total_seconds() > 120: # 2분으로 상향
                logger.warning(f"🧟 Zombie Job Detected: {job_id}. Marking as failed.")
                await update_job_status(job_id, status="failed", errors=["Zombie process detected (timeout)"])
                job["status"] = "failed"
        except Exception as e:
            logger.error(f"Error parsing started_at: {e}")
            
    return job

@router.post("/sync/{job_id}/stop")
async def stop_sync_job(job_id: str):
    await update_job_status(job_id, cancel_requested=True)
    return {"message": "Stop request sent to DB."}

@router.post("/workers/reset")
async def reset_all_workers():
    """모든 워커 강제 중단 및 잠금 해제"""
    try:
        # 1. 모든 잠금 해제
        supabase.table("job_locks").delete().neq("id", "00000000-0000-0000-0000-000000000000").execute()
        
        # 2. 실행 중인 작업 상태 변경
        # 2. 실행 중인 작업 상태 변경 (먼저 취소 요청 플래그 설정)
        supabase.table("sync_jobs").update({
            "cancel_requested": True
        }).eq("status", "running").execute()

        # 3. 잠시 대기 후 강제 종료 처리 (옵션: 워커가 스스로 종료하도록 유도하거나, 여기서 강제로 status를 바꿈)
        # 여기서는 즉시 status를 stopped로 바꾸지만, 워커 로직에서 cancel_requested를 확인하고 종료하도록 함
        supabase.table("sync_jobs").update({
            "status": "stopped",
            "errors": ["Force stopped by admin"],
            "completed_at": datetime.utcnow().isoformat()
        }).eq("status", "running").execute()
        
        return {"message": "All workers have been reset and locks released."}
    except Exception as e:
        logger.error(f"Failed to reset workers: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# --- Bulk Import & AI Refiner Endpoints ---

@router.post("/bulk/import", response_model=SyncJobResponse)
async def run_bulk_import(background_tasks: BackgroundTasks, max_studies: int = 5000, mode: str = "daily"):
    """
    ClinicalTrials.gov 전체 덤프에서 ADC 데이터 일괄 임포트
    mode: 'daily' (기본값) 또는 'full' (전체 적재)
    """
    from app.services.bulk_importer import BulkImporter
    
    # 중복 실행 방지
    if not await job_lock.acquire("bulk_import"):
        raise HTTPException(
            status_code=409, 
            detail="이미 실행 중인 Bulk Import 작업이 있습니다. 나중에 다시 시도하세요."
        )
    
    job_id = f"bulk_import_{uuid4().hex[:8]}"
    data = {
        "id": job_id, 
        "status": "queued", 
        "source": "clinical_trials_bulk", 
        "started_at": datetime.utcnow().isoformat()
    }
    supabase.table("sync_jobs").insert(data).execute()
    
    async def run_with_lock_release(importer, job_id, max_studies, mode):
        try:
            await importer.run_import(job_id, max_studies, mode=mode)
        finally:
            await job_lock.release("bulk_import")
    
    importer = BulkImporter()
    background_tasks.add_task(run_with_lock_release, importer, job_id, max_studies, mode)
    
    return SyncJobResponse(
        job_id=job_id, 
        status="queued", 
        message=f"Bulk import started (mode: {mode}, max {max_studies} studies)."
    )

@router.post("/refiner/run", response_model=SyncJobResponse)
async def run_ai_refiner(background_tasks: BackgroundTasks, max_records: int = 50):
    """미정제 레코드 LLM 분석 및 SMILES 보강"""
    from app.services.ai_refiner import ai_refiner
    
    # 중복 실행 방지
    if not await job_lock.acquire("ai_refiner"):
        raise HTTPException(
            status_code=409, 
            detail="이미 실행 중인 AI Refiner 작업이 있습니다. 나중에 다시 시도하세요."
        )
    
    job_id = f"ai_refiner_{uuid4().hex[:8]}"
    data = {
        "id": job_id, 
        "status": "queued", 
        "source": "ai_refiner", 
        "started_at": datetime.utcnow().isoformat()
    }
    supabase.table("sync_jobs").insert(data).execute()
    
    async def run_with_lock_release(job_id, max_records):
        try:
            await ai_refiner.process_pending_records(job_id, max_records)
        finally:
            await job_lock.release("ai_refiner")
    
    background_tasks.add_task(run_with_lock_release, job_id, max_records)
    
    return SyncJobResponse(
        job_id=job_id, 
        status="queued", 
        message=f"AI Refiner started (max {max_records} records)."
    )

@router.post("/sync/pubmed-knowledge", response_model=SyncJobResponse)
async def sync_pubmed_knowledge(
    background_tasks: BackgroundTasks, 
    batch_size: int = 100,
    mode: str = "incremental"
):
    """
    PubMed 기반 ADC Knowledge Base 구축
    golden_set_library에서 약물 추출 → PubMed 검색 → Gemini 분석 → knowledge_base 저장
    
    mode: 'incremental' (배치 크기만큼) 또는 'full' (전체)
    """
    from app.services.pubmed_knowledge_service import PubMedKnowledgeService
    
    # 중복 실행 방지
    if not await job_lock.acquire("pubmed_knowledge"):
        raise HTTPException(
            status_code=409, 
            detail="이미 실행 중인 PubMed Knowledge 작업이 있습니다."
        )
    
    job_id = f"pubmed_knowledge_{uuid4().hex[:8]}"
    data = {
        "id": job_id, 
        "status": "queued", 
        "source": "pubmed_knowledge", 
        "started_at": datetime.utcnow().isoformat()
    }
    supabase.table("sync_jobs").insert(data).execute()
    
    async def run_with_lock_release(job_id, batch_size, mode):
        try:
            service = PubMedKnowledgeService()
            await service.run_batch(job_id, batch_size, mode)
        finally:
            await job_lock.release("pubmed_knowledge")
    
    background_tasks.add_task(run_with_lock_release, job_id, batch_size, mode)
    
    return SyncJobResponse(
        job_id=job_id, 
        status="queued", 
        message=f"PubMed Knowledge sync started (batch: {batch_size}, mode: {mode})."
    )

@router.get("/refiner/pending")
async def get_pending_count():
    """미정제 레코드 수 조회"""
    try:
        res = supabase.table("golden_set_library")\
            .select("count", count="exact")\
            .eq("ai_refined", False)\
            .execute()
        return {"pending_count": res.count}
    except Exception as e:
        return {"pending_count": 0, "error": str(e)}

@router.get("/health")
async def health():
    res = supabase.table("sync_jobs").select("count", count="exact").eq("status", "running").execute()
    return {"status": "healthy", "active_jobs": res.count}
