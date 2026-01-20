"""
Data Worker API - 외부 데이터 소스 동기화
ClinicalTrials.gov, PubMed 등에서 데이터 수집 후 draft 상태로 저장
"""
from fastapi import APIRouter, BackgroundTasks, HTTPException
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
from datetime import datetime, timedelta
from uuid import uuid4
import httpx
import asyncio
from Bio import Entrez
from tenacity import retry, stop_after_attempt, wait_exponential
import logging
import json
import re

from app.services.perplexity_service import enrich_drug_data
from app.core.config import settings
from app.core.supabase import supabase
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from app.services.openfda_service import openfda_service

router = APIRouter()
logger = logging.getLogger(__name__)

# Entrez 설정
Entrez.email = settings.NCBI_EMAIL
if settings.NCBI_API_KEY:
    Entrez.api_key = settings.NCBI_API_KEY
if settings.NCBI_TOOL:
    Entrez.tool = settings.NCBI_TOOL

# In-memory 동기화 상태 추적 (TODO: Redis로 대체)
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
    """
    LLM을 사용하여 약물 이름, 타겟, 페이로드, 링커, 독성 정보 및 근거 추출
    """
    try:
        llm_refiner = ChatOpenAI(
            model=settings.FAST_LLM, # gpt-4o-mini
            temperature=0,
            api_key=settings.OPENAI_API_KEY
        )
        
        system_prompt = """You are an expert Clinical & Pharmaceutical Data Analyst. 
        Extract structured ADC (Antibody-Drug Conjugate) information from the provided text.
        
        For each field, provide the 'value' and the 'evidence' (the exact sentence or phrase from the text).
        
        Fields to extract:
        1. **drug_name**: Specific ADC code name (e.g., DS-8201, IMGN-632) or generic name.
        2. **target**: The biological target (e.g., HER2, TROP2, CD33).
        3. **payload**: The cytotoxic drug (e.g., MMAE, DXd, SN-38).
        4. **linker**: The linker type (e.g., Val-Cit, tetrapeptide-based).
        5. **outcome_type**: 'Success', 'Failure', or 'Ongoing'.
        6. **failure_reason**: If 'Failure', summarize the scientific reason (e.g., "Due to Grade 3 ILD").
        7. **toxicity_profile**: Key safety issues or boxed warnings.
        
        Output JSON format:
        {
            "drug_name": {"value": "DS-8201a", "evidence": "...fam-trastuzumab deruxtecan-nxki..."},
            "target": {"value": "HER2", "evidence": "...targeting human epidermal growth factor receptor 2..."},
            "payload": {"value": "DXd", "evidence": "...deruxtecan, a derivative of exatecan..."},
            "linker": {"value": "Tetrapeptide-based", "evidence": "...cleavable tetrapeptide-based linker..."},
            "outcome_type": {"value": "Success", "evidence": "...demonstrated significant improvement in PFS..."},
            "failure_reason": {"value": null, "evidence": null},
            "toxicity_profile": {"value": "ILD/Pneumonitis", "evidence": "...Boxed Warning for interstitial lung disease..."}
        }
        """
        
        user_prompt = f"""
        Title: {title}
        Content: {text_content}
        Status Context: {raw_status if raw_status else 'N/A'}
        """
        
        prompt = ChatPromptTemplate.from_messages([
            ("system", system_prompt),
            ("user", user_prompt)
        ])
    
        chain = prompt | llm_refiner
        response = await chain.ainvoke({})
        
        content = response.content.strip()
        if "```json" in content:
            content = content.split("```json")[1].split("```")[0]
        elif "```" in content:
            content = content.split("```")[1].split("```")[0]
            
        return json.loads(content)
        
    except Exception as e:
        logger.error(f"LLM Refine Error: {e}")
        return {
            "drug_name": {"value": "Unknown", "evidence": None},
            "target": {"value": "Unknown", "evidence": None},
            "payload": {"value": "Unknown", "evidence": None},
            "linker": {"value": "Unknown", "evidence": None},
            "outcome_type": {"value": "Unknown", "evidence": None},
            "failure_reason": {"value": None, "evidence": None},
            "toxicity_profile": {"value": "Unknown", "evidence": None}
        }

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
async def fetch_clinical_trials(query: str = "ADC OR Antibody-Drug Conjugate"):
    """ClinicalTrials.gov API v2 호출"""
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.term": query,
        "pageSize": 50,
        "format": "json"
    }
    async with httpx.AsyncClient() as client:
        response = await client.get(url, params=params, timeout=30.0)
        response.raise_for_status()
        return response.json()

async def process_clinical_trials_data(job_id: str):
    """ClinicalTrials.gov 데이터 수집 및 처리"""
    if job_id not in sync_jobs:
        sync_jobs[job_id] = {"status": "queued", "records_found": 0, "records_drafted": 0, "errors": []}
    
    sync_jobs[job_id]["status"] = "running"
    
    try:
        data = await fetch_clinical_trials()
        studies = data.get("studies", [])
        sync_jobs[job_id]["records_found"] = len(studies)
        drafted = 0
        
        for study in studies:
            try:
                # 중단 요청 체크
                if sync_jobs.get(job_id, {}).get("cancel_requested"):
                    sync_jobs[job_id]["status"] = "stopped"
                    return

                protocol = study.get("protocolSection", {})
                id_module = protocol.get("identificationModule", {})
                nct_id = id_module.get("nctId")
                
                existing = supabase.table("golden_set_library").select("id").eq("properties->>nct_id", nct_id).execute()
                if existing.data:
                    continue

                conditions = protocol.get("conditionsModule", {}).get("conditions", [])
                interventions = protocol.get("armsInterventionsModule", {}).get("interventions", [])
                
                drug_name = "Unknown"
                for intervention in interventions:
                    if intervention.get("type") == "DRUG":
                        drug_name = intervention.get("name")
                        break
                
                refined_data = await refine_drug_data_with_llm(
                    title=id_module.get("officialTitle"),
                    text_content=str(protocol.get("descriptionModule", {})),
                    raw_status=protocol.get("statusModule", {}).get("overallStatus")
                )
                
                final_name = refined_data.get("drug_name", {}).get("value", drug_name)
                if final_name == "Unknown": final_name = drug_name
                
                outcome_type = refined_data.get("outcome_type", {}).get("value", "Unknown")
                overall_status = protocol.get("statusModule", {}).get("overallStatus", "Unknown")
                
                if outcome_type == "Unknown":
                    if overall_status in ["COMPLETED", "APPROVED"]: outcome_type = "Success"
                    elif overall_status in ["TERMINATED", "WITHDRAWN", "SUSPENDED"]: outcome_type = "Failure"
                    elif overall_status in ["RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION"]: outcome_type = "Ongoing"

                new_entry = {
                    "name": final_name,
                    "category": "clinical_trial",
                    "description": f"[{outcome_type}] {id_module.get('officialTitle')}",
                    "properties": {
                        "nct_id": nct_id,
                        "title": id_module.get("officialTitle"),
                        "conditions": conditions,
                        "drug_name": drug_name,
                        "phase": protocol.get("designModule", {}).get("phases", ["Unknown"])[0],
                        "ai_analysis": refined_data, 
                        "target": refined_data.get("target", {}).get("value"),
                        "payload": refined_data.get("payload", {}).get("value"),
                        "linker": refined_data.get("linker", {}).get("value"),
                        "toxicity_profile": refined_data.get("toxicity_profile", {}).get("value"),
                        "overall_status": overall_status,
                        "why_stopped": protocol.get("statusModule", {}).get("whyStopped")
                    },
                    "status": "draft",
                    "enrichment_source": "clinical_trials_ai_refined",
                    "raw_data": study,
                    "outcome_type": outcome_type,
                    "failure_reason": refined_data.get("failure_reason", {}).get("value") if outcome_type == "Failure" else None,
                    "is_ai_extracted": refined_data.get("drug_name", {}).get("value") != "Unknown"
                }

                supabase.table("golden_set_library").insert(new_entry).execute()
                drafted += 1
                sync_jobs[job_id]["records_drafted"] = drafted
                
            except Exception as e:
                logger.error(f"Clinical Trial Error ({nct_id}): {e}")
                sync_jobs[job_id]["errors"].append(str(e))
        
        sync_jobs[job_id]["status"] = "completed"
        sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
                    
    except Exception as e:
        logger.error(f"Clinical Trial Worker Error: {e}")
        sync_jobs[job_id]["status"] = "failed"
        sync_jobs[job_id]["errors"].append(str(e))

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def fetch_pubmed_ids(term: str = "Antibody-Drug Conjugate", max_results: int = 50):
    """PubMed ID 검색 (최신순 정렬)"""
    handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results, sort="pub_date")
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def fetch_pubmed_details(id_list: List[str]):
    """PubMed 상세 정보 조회"""
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    return records

async def process_pubmed_data(job_id: str):
    """PubMed 데이터 수집"""
    if job_id not in sync_jobs:
        sync_jobs[job_id] = {"status": "queued", "records_found": 0, "records_drafted": 0, "errors": []}
    sync_jobs[job_id]["status"] = "running"
    
    try:
        loop = asyncio.get_event_loop()
        search_term = "Antibody-Drug Conjugate OR ADC"
        id_list = await loop.run_in_executor(None, fetch_pubmed_ids, search_term)
        
        sync_jobs[job_id]["records_found"] = len(id_list)
        if not id_list:
            sync_jobs[job_id]["status"] = "completed"
            sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
            return

        papers = await loop.run_in_executor(None, fetch_pubmed_details, id_list)
        drafted = 0
        
        if 'PubmedArticle' in papers:
            for article in papers['PubmedArticle']:
                try:
                    # 중단 요청 체크
                    if sync_jobs.get(job_id, {}).get("cancel_requested"):
                        sync_jobs[job_id]["status"] = "stopped"
                        return

                    medline = article['MedlineCitation']
                    article_data = medline['Article']
                    pmid = str(medline['PMID'])
                    title = article_data.get('ArticleTitle', 'No Title')
                    
                    abstract_texts = article_data.get('Abstract', {}).get('AbstractText', [])
                    if isinstance(abstract_texts, list):
                        abstract = " ".join([str(t) for t in abstract_texts])
                    else:
                        abstract = str(abstract_texts)
                    
                    if not abstract or abstract == "None":
                        abstract = "No Abstract available."

                    existing = supabase.table("knowledge_base").select("id").eq("title", title).execute()
                    if existing.data:
                        continue
                        
                    refined_data = await refine_drug_data_with_llm(title=title, text_content=abstract)
                    new_kb = {
                        "source_type": "PubMed",
                        "title": title,
                        "summary": abstract[:500] + "...",
                        "content": abstract,
                        "relevance_score": 0.0,
                        "source_tier": 1,
                        "rag_status": "pending",
                        "ai_analysis": refined_data
                    }
                    supabase.table("knowledge_base").insert(new_kb).execute()
                    drafted += 1
                    sync_jobs[job_id]["records_drafted"] = drafted
                except Exception as e:
                    logger.error(f"PubMed Article Error (PMID {pmid}): {e}")
                    sync_jobs[job_id]["errors"].append(f"PMID {pmid}: {str(e)}")

        sync_jobs[job_id]["status"] = "completed"
        sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    except Exception as e:
        logger.error(f"PubMed Worker Error: {e}")
        sync_jobs[job_id]["status"] = "failed"
        sync_jobs[job_id]["errors"].append(str(e))

async def wrapper_openfda_sync(job_id: str):
    """OpenFDA 서비스 실행 및 상태 업데이트 래퍼"""
    sync_jobs[job_id]["status"] = "running"
    try:
        await openfda_service.sync_to_db(job_id=job_id)
        if sync_jobs[job_id]["status"] != "stopped":
            sync_jobs[job_id]["status"] = "completed"
            sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    except Exception as e:
        logger.error(f"OpenFDA Wrapper Error: {e}")
        sync_jobs[job_id]["status"] = "failed"
        sync_jobs[job_id]["errors"].append(str(e))

async def wrapper_creative_crawl(job_id: str, search_term: Optional[str] = None):
    """Creative Biolabs 크롤러 실행 및 상태 업데이트 래퍼"""
    from app.services.creative_biolabs_crawler import creative_crawler
    sync_jobs[job_id]["status"] = "running"
    try:
        await creative_crawler.run(search_term=search_term, job_id=job_id)
        if sync_jobs[job_id]["status"] != "stopped":
            sync_jobs[job_id]["status"] = "completed"
            sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
    except Exception as e:
        logger.error(f"Creative Crawler Wrapper Error: {e}")
        sync_jobs[job_id]["status"] = "failed"
        sync_jobs[job_id]["errors"].append(str(e))

@router.get("/settings/{source_id}")
async def get_source_settings(source_id: str):
    try:
        res = supabase.table("data_source_settings").select("*").eq("source_id", source_id).execute()
        if res.data:
            return res.data[0]
        return {"source_id": source_id, "auto_sync": False, "sync_interval_hours": 24}
    except Exception:
        return {"source_id": source_id, "auto_sync": False, "sync_interval_hours": 24}

@router.post("/settings")
async def save_source_settings(setting: DataSourceSetting):
    try:
        data = {
            "source_id": setting.source_id,
            "auto_sync": setting.auto_sync,
            "sync_interval_hours": setting.sync_interval_hours,
            "updated_at": datetime.utcnow().isoformat()
        }
        supabase.table("data_source_settings").upsert(data).execute()
        return {"status": "success", "message": f"Settings for {setting.source_id} saved."}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/sync/clinical", response_model=SyncJobResponse)
async def sync_clinical_trials(background_tasks: BackgroundTasks):
    job_id = f"sync_clinical_{uuid4().hex[:8]}"
    sync_jobs[job_id] = {
        "status": "queued", "source": "clinical_trials", "records_found": 0, "records_drafted": 0,
        "started_at": datetime.utcnow().isoformat(), "completed_at": None, "cancel_requested": False, "errors": []
    }
    background_tasks.add_task(process_clinical_trials_data, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="ClinicalTrials.gov sync started.")

@router.post("/sync/pubmed", response_model=SyncJobResponse)
async def sync_pubmed(background_tasks: BackgroundTasks):
    job_id = f"sync_pubmed_{uuid4().hex[:8]}"
    sync_jobs[job_id] = {
        "status": "queued", "source": "pubmed", "records_found": 0, "records_drafted": 0,
        "started_at": datetime.utcnow().isoformat(), "completed_at": None, "cancel_requested": False, "errors": []
    }
    background_tasks.add_task(process_pubmed_data, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="PubMed sync started.")

@router.post("/sync/openfda", response_model=SyncJobResponse)
async def sync_openfda(background_tasks: BackgroundTasks):
    job_id = f"sync_openfda_{uuid4().hex[:8]}"
    sync_jobs[job_id] = {
        "status": "queued", "source": "openfda", "records_found": 0, "records_drafted": 0,
        "started_at": datetime.utcnow().isoformat(), "completed_at": None, "cancel_requested": False, "errors": []
    }
    background_tasks.add_task(wrapper_openfda_sync, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="OpenFDA Approved ADC sync started.")

@router.post("/crawler/creative/run", response_model=SyncJobResponse)
async def run_creative_crawler(background_tasks: BackgroundTasks, search_term: Optional[str] = None):
    job_id = f"crawl_creative_{uuid4().hex[:8]}"
    sync_jobs[job_id] = {
        "status": "queued", "source": "creative_biolabs", "records_found": 0, "records_drafted": 0,
        "started_at": datetime.utcnow().isoformat(), "completed_at": None, "cancel_requested": False, "errors": []
    }
    background_tasks.add_task(wrapper_creative_crawl, job_id, search_term)
    return SyncJobResponse(job_id=job_id, status="queued", message="Creative Biolabs crawler started.")

@router.post("/sync/goldenset", response_model=SyncJobResponse)
async def sync_goldenset(background_tasks: BackgroundTasks):
    job_id = f"refresh_{uuid4().hex[:8]}"
    sync_jobs[job_id] = {
        "status": "completed", "source": "goldenset", "records_found": 0, "records_drafted": 0,
        "started_at": datetime.utcnow().isoformat(), "completed_at": datetime.utcnow().isoformat(),
        "cancel_requested": False, "errors": []
    }
    return SyncJobResponse(job_id=job_id, status="completed", message="Golden Set library view refreshed.")

@router.get("/sync/{job_id}", response_model=SyncJobStatus)
async def get_sync_status(job_id: str):
    if job_id not in sync_jobs:
        raise HTTPException(status_code=404, detail=f"Sync job {job_id} not found")
    return SyncJobStatus(job_id=job_id, **sync_jobs[job_id])

@router.post("/sync/{job_id}/stop")
async def stop_sync_job(job_id: str):
    if job_id not in sync_jobs:
        raise HTTPException(status_code=404, detail="Sync job not found")
    sync_jobs[job_id]["cancel_requested"] = True
    return {"message": "Stop request sent to worker."}

@router.get("/health")
async def scheduler_health():
    return {
        "status": "healthy",
        "active_jobs": len([j for j in sync_jobs.values() if j.get("status") == "running"])
    }
