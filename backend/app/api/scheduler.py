"""
Data Worker API - 외부 데이터 소스 동기화 (AstraForge 2.0)
비동기 병렬 처리(Asyncio) 및 자동 스케줄링(APScheduler) 지원
"""
from fastapi import APIRouter, BackgroundTasks, HTTPException
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
from datetime import datetime
from uuid import uuid4
import httpx
import asyncio
from Bio import Entrez
from tenacity import retry, stop_after_attempt, wait_exponential
import logging
import json
import re

from app.core.config import settings
from app.core.supabase import supabase
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from app.services.openfda_service import openfda_service
from app.services.scheduler_engine import scheduler_engine

router = APIRouter()
logger = logging.getLogger(__name__)

# Entrez 설정
Entrez.email = settings.NCBI_EMAIL
if settings.NCBI_API_KEY:
    Entrez.api_key = settings.NCBI_API_KEY
if settings.NCBI_TOOL:
    Entrez.tool = settings.NCBI_TOOL

# In-memory 동기화 상태 추적
sync_jobs: Dict[str, Dict[str, Any]] = {}

class DataSourceSetting(BaseModel):
    source_id: str
    auto_sync: bool
    sync_interval_hours: int

class SyncJobResponse(BaseModel):
    job_id: str
    status: str
    message: str

class SyncJobStatus(BaseModel):
    job_id: str
    status: str
    source: str
    records_found: int
    records_drafted: int
    started_at: str
    completed_at: Optional[str] = None
    errors: List[str] = []

async def refine_drug_data_with_llm(title: str, text_content: str, raw_status: Optional[str] = None) -> dict:
    """LLM을 사용하여 ADC 정보 추출"""
    try:
        llm_refiner = ChatOpenAI(
            model=settings.FAST_LLM,
            temperature=0,
            api_key=settings.OPENAI_API_KEY
        )
        
        system_prompt = """You are an expert Clinical & Pharmaceutical Data Analyst. 
        Extract structured ADC (Antibody-Drug Conjugate) information from the provided text.
        Output JSON format:
        {
            "drug_name": {"value": "DS-8201a", "evidence": "..."},
            "target": {"value": "HER2", "evidence": "..."},
            "payload": {"value": "DXd", "evidence": "..."},
            "linker": {"value": "Tetrapeptide-based", "evidence": "..."},
            "outcome_type": {"value": "Success", "evidence": "..."},
            "failure_reason": {"value": null, "evidence": null},
            "toxicity_profile": {"value": "ILD/Pneumonitis", "evidence": "..."}
        }"""
        
        user_prompt = f"Title: {title}\nContent: {text_content}\nStatus: {raw_status}"
        prompt = ChatPromptTemplate.from_messages([("system", system_prompt), ("user", user_prompt)])
        chain = prompt | llm_refiner
        response = await chain.ainvoke({})
        content = response.content.strip()
        if "```json" in content: content = content.split("```json")[1].split("```")[0]
        return json.loads(content)
    except Exception as e:
        logger.error(f"LLM Refine Error: {e}")
        return {"drug_name": {"value": "Unknown"}}

# --- ClinicalTrials.gov Worker (Async & Paging) ---

async def process_clinical_trials_data(job_id: str, max_records: int = 1000):
    """ClinicalTrials.gov 대량 수집 워커"""
    if job_id not in sync_jobs:
        sync_jobs[job_id] = {"status": "running", "source": "clinical_trials", "records_found": 0, "records_drafted": 0, "errors": [], "started_at": datetime.utcnow().isoformat()}
    
    sync_jobs[job_id]["status"] = "running"
    base_url = "https://clinicaltrials.gov/api/v2/studies"
    next_token = None
    total_collected = 0
    drafted = 0

    async with httpx.AsyncClient(timeout=30.0) as client:
        while total_collected < max_records:
            if sync_jobs.get(job_id, {}).get("cancel_requested"):
                sync_jobs[job_id]["status"] = "stopped"
                return

            params = {"query.term": "ADC OR Antibody-Drug Conjugate", "pageSize": 50, "format": "json"}
            if next_token: params["pageToken"] = next_token

            try:
                res = await client.get(base_url, params=params)
                res.raise_for_status()
                data = res.json()
                studies = data.get("studies", [])
                if not studies: break

                sync_jobs[job_id]["records_found"] += len(studies)
                
                # 병렬 처리: 5개씩 묶어서 처리 (LLM API Rate Limit 고려)
                for i in range(0, len(studies), 5):
                    batch = studies[i:i+5]
                    tasks = [process_single_study(study, job_id) for study in batch]
                    results = await asyncio.gather(*tasks)
                    drafted += sum(1 for r in results if r)
                    sync_jobs[job_id]["records_drafted"] = drafted
                    
                    if sync_jobs.get(job_id, {}).get("cancel_requested"):
                        sync_jobs[job_id]["status"] = "stopped"
                        return

                total_collected += len(studies)
                next_token = data.get("nextPageToken")
                if not next_token: break
                
                await asyncio.sleep(1) # 매너 딜레이
            except Exception as e:
                logger.error(f"ClinicalTrials Batch Error: {e}")
                sync_jobs[job_id]["errors"].append(str(e))
                break

    sync_jobs[job_id]["status"] = "completed"
    sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()

async def process_single_study(study: dict, job_id: str) -> bool:
    """개별 스터디 처리 및 저장"""
    try:
        protocol = study.get("protocolSection", {})
        id_module = protocol.get("identificationModule", {})
        nct_id = id_module.get("nctId")
        
        # 중복 체크
        existing = supabase.table("golden_set_library").select("id").eq("properties->>nct_id", nct_id).execute()
        if existing.data: return False

        refined = await refine_drug_data_with_llm(
            title=id_module.get("officialTitle"),
            text_content=str(protocol.get("descriptionModule", {})),
            raw_status=protocol.get("statusModule", {}).get("overallStatus")
        )
        
        outcome_type = refined.get("outcome_type", {}).get("value", "Unknown")
        new_entry = {
            "name": refined.get("drug_name", {}).get("value", "Unknown"),
            "category": "clinical_trial",
            "description": id_module.get("officialTitle"),
            "properties": {"nct_id": nct_id, "ai_analysis": refined},
            "status": "draft",
            "outcome_type": outcome_type,
            "enrichment_source": "clinical_trials_v2_async"
        }
        supabase.table("golden_set_library").insert(new_entry).execute()
        return True
    except Exception as e:
        logger.error(f"Study Error {nct_id}: {e}")
        return False

# --- PubMed Worker (Async & Paging) ---

async def process_pubmed_data(job_id: str, max_records: int = 500):
    """PubMed 대량 수집 워커"""
    if job_id not in sync_jobs:
        sync_jobs[job_id] = {"status": "running", "source": "pubmed", "records_found": 0, "records_drafted": 0, "errors": [], "started_at": datetime.utcnow().isoformat()}
    
    sync_jobs[job_id]["status"] = "running"
    drafted = 0
    
    try:
        loop = asyncio.get_event_loop()
        search_term = "Antibody-Drug Conjugate OR ADC OR 'ADC drug' OR 'ADC therapy'"
        id_list = await loop.run_in_executor(None, lambda: fetch_pubmed_ids(search_term, max_records))
        sync_jobs[job_id]["records_found"] = len(id_list)
        
        # 20개씩 끊어서 병렬 상세 조회 및 처리
        for i in range(0, len(id_list), 20):
            if sync_jobs.get(job_id, {}).get("cancel_requested"):
                sync_jobs[job_id]["status"] = "stopped"
                return

            batch_ids = id_list[i:i+20]
            papers = await loop.run_in_executor(None, fetch_pubmed_details, batch_ids)
            
            if 'PubmedArticle' in papers:
                tasks = [process_single_article(article, job_id) for article in papers['PubmedArticle']]
                results = await asyncio.gather(*tasks)
                drafted += sum(1 for r in results if r)
                sync_jobs[job_id]["records_drafted"] = drafted

        sync_jobs[job_id]["status"] = "completed"
        sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    except Exception as e:
        logger.error(f"PubMed Worker Error: {e}")
        sync_jobs[job_id]["status"] = "failed"
        sync_jobs[job_id]["errors"].append(str(e))

async def process_single_article(article: dict, job_id: str) -> bool:
    """개별 논문 처리 및 저장"""
    try:
        medline = article['MedlineCitation']
        article_data = medline['Article']
        title = article_data.get('ArticleTitle', 'No Title')
        
        existing = supabase.table("knowledge_base").select("id").eq("title", title).execute()
        if existing.data: return False

        abstract_texts = article_data.get('Abstract', {}).get('AbstractText', [])
        abstract = " ".join([str(t) for t in abstract_texts]) if isinstance(abstract_texts, list) else str(abstract_texts)
        
        refined = await refine_drug_data_with_llm(title=title, text_content=abstract)
        new_kb = {
            "source_type": "PubMed",
            "title": title,
            "content": abstract,
            "ai_analysis": refined,
            "rag_status": "pending"
        }
        supabase.table("knowledge_base").insert(new_kb).execute()
        return True
    except Exception:
        return False

def fetch_pubmed_ids(term: str, max_results: int):
    handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results, sort="pub_date")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_pubmed_details(id_list: List[str]):
    handle = Entrez.efetch(db="pubmed", id=",".join(id_list), retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

# --- API Endpoints ---

@router.post("/settings")
async def save_source_settings(setting: DataSourceSetting):
    try:
        data = {"source_id": setting.source_id, "auto_sync": setting.auto_sync, "sync_interval_hours": setting.sync_interval_hours, "updated_at": datetime.utcnow().isoformat()}
        supabase.table("data_source_settings").upsert(data).execute()
        
        # 스케줄러 업데이트
        if setting.auto_sync:
            job_func = None
            if setting.source_id == "clinical_trials": job_func = process_clinical_trials_data
            elif setting.source_id == "pubmed": job_func = process_pubmed_data
            
            if job_func:
                scheduler_engine.add_or_update_job(f"auto_{setting.source_id}", job_func, setting.sync_interval_hours, [f"auto_{setting.source_id}"])
        else:
            scheduler_engine.remove_job(f"auto_{setting.source_id}")
            
        return {"status": "success"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/sync/clinical", response_model=SyncJobResponse)
async def sync_clinical_trials(background_tasks: BackgroundTasks):
    job_id = f"sync_clinical_{uuid4().hex[:8]}"
    sync_jobs[job_id] = {"status": "queued", "source": "clinical_trials", "records_found": 0, "records_drafted": 0, "started_at": datetime.utcnow().isoformat(), "cancel_requested": False, "errors": []}
    background_tasks.add_task(process_clinical_trials_data, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="ClinicalTrials.gov sync started.")

@router.post("/sync/pubmed", response_model=SyncJobResponse)
async def sync_pubmed(background_tasks: BackgroundTasks):
    job_id = f"sync_pubmed_{uuid4().hex[:8]}"
    sync_jobs[job_id] = {"status": "queued", "source": "pubmed", "records_found": 0, "records_drafted": 0, "started_at": datetime.utcnow().isoformat(), "cancel_requested": False, "errors": []}
    background_tasks.add_task(process_pubmed_data, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="PubMed sync started.")

@router.post("/sync/openfda", response_model=SyncJobResponse)
async def sync_openfda(background_tasks: BackgroundTasks):
    job_id = f"sync_openfda_{uuid4().hex[:8]}"
    sync_jobs[job_id] = {"status": "queued", "source": "openfda", "records_found": 0, "records_drafted": 0, "started_at": datetime.utcnow().isoformat(), "cancel_requested": False, "errors": []}
    background_tasks.add_task(openfda_service.sync_to_db, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="OpenFDA sync started.")

@router.post("/crawler/creative/run", response_model=SyncJobResponse)
async def run_creative_crawler(background_tasks: BackgroundTasks, search_term: Optional[str] = None):
    from app.services.creative_biolabs_crawler import creative_crawler
    job_id = f"crawl_creative_{uuid4().hex[:8]}"
    sync_jobs[job_id] = {"status": "queued", "source": "creative_biolabs", "records_found": 0, "records_drafted": 0, "started_at": datetime.utcnow().isoformat(), "cancel_requested": False, "errors": []}
    background_tasks.add_task(creative_crawler.run, search_term, 3, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="Creative Biolabs crawler started.")

@router.get("/sync/{job_id}", response_model=SyncJobStatus)
async def get_sync_status(job_id: str):
    if job_id not in sync_jobs: raise HTTPException(status_code=404, detail="Not found")
    return SyncJobStatus(job_id=job_id, **sync_jobs[job_id])

@router.post("/sync/{job_id}/stop")
async def stop_sync_job(job_id: str):
    if job_id in sync_jobs: sync_jobs[job_id]["cancel_requested"] = True
    return {"message": "Stop request sent."}

@router.get("/health")
async def health():
    return {"status": "healthy", "active_jobs": len([j for j in sync_jobs.values() if j.get("status") == "running"])}
