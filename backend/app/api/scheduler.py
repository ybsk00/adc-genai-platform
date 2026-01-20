"""
Data Worker API - 외부 데이터 소스 동기화 (AstraForge 2.0)
DB 기반 상태 관리(Supabase) 및 PubChem 화학 정보 연동
"""
from fastapi import APIRouter, BackgroundTasks, HTTPException
from pydantic import BaseModel
from typing import Optional, List, Dict, Any
from datetime import datetime
from uuid import uuid4
import aiohttp
import asyncio
from Bio import Entrez
import logging
import json

from app.core.config import settings
from app.core.supabase import supabase
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
from app.services.openfda_service import openfda_service
from app.services.scheduler_engine import scheduler_engine
from app.services.chemical_resolver import chemical_resolver

router = APIRouter()
logger = logging.getLogger(__name__)

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
        
        supabase.table("sync_jobs").update(kwargs).eq("id", job_id).execute()
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
    """ClinicalTrials.gov 대량 수집 워커 (DB 기반)"""
    await update_job_status(job_id, status="running")
    
    base_url = "https://clinicaltrials.gov/api/v2/studies"
    next_token = None
    total_collected = 0
    drafted = 0
    errors = []
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
        'Accept': 'application/json'
    }

    try:
        async with aiohttp.ClientSession(headers=headers) as session:
            while total_collected < max_records:
                if await is_cancelled(job_id):
                    await update_job_status(job_id, status="stopped")
                    return

                params = {"query.term": "antibody drug conjugate", "pageSize": 50, "format": "json"}
                if next_token: params["pageToken"] = next_token

                try:
                    async with session.get(base_url, params=params, timeout=aiohttp.ClientTimeout(total=60)) as res:
                        if res.status == 403:
                            errors.append("ClinicalTrials API returned 403 Forbidden - API access restricted")
                            logger.error("ClinicalTrials API 403 Forbidden")
                            break
                        if res.status != 200:
                            errors.append(f"ClinicalTrials API error: {res.status}")
                            break
                        
                        data = await res.json()
                        studies = data.get("studies", [])
                        if not studies: break

                        # 발견 즉시 DB 업데이트
                        job = await get_job_from_db(job_id)
                        current_found = job.get("records_found", 0) if job else 0
                        await update_job_status(job_id, records_found=current_found + len(studies))
                        
                        # 단순 저장 (LLM 없이)
                        for study in studies:
                            if await is_cancelled(job_id):
                                await update_job_status(job_id, status="stopped")
                                return
                            try:
                                result = await process_single_study_simple(study)
                                if result: drafted += 1
                            except Exception as e:
                                errors.append(str(e)[:100])
                        
                        await update_job_status(job_id, records_drafted=drafted, errors=errors[:20])

                        total_collected += len(studies)
                        next_token = data.get("nextPageToken")
                        if not next_token: break
                        await asyncio.sleep(0.5)
                except Exception as e:
                    errors.append(f"Request error: {str(e)[:100]}")
                    logger.error(f"ClinicalTrials request error: {e}")
                    break

        await update_job_status(job_id, status="completed", completed_at=datetime.utcnow().isoformat(), errors=errors[:20])
    except Exception as e:
        logger.error(f"ClinicalTrials Global Error: {e}")
        await update_job_status(job_id, status="failed", errors=[str(e)])

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

async def process_pubmed_data(job_id: str, max_records: int = 500):
    """PubMed 대량 수집 워커 (DB 기반)"""
    await update_job_status(job_id, status="running")
    drafted = 0
    errors = []
    
    try:
        loop = asyncio.get_event_loop()
        search_term = "Antibody-Drug Conjugate OR ADC OR ADC drug OR ADC therapy"
        
        id_list = await loop.run_in_executor(None, lambda: fetch_pubmed_ids(search_term, max_records))
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

def fetch_pubmed_ids(term: str, max_results: int):
    handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results, sort="pub_date")
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
async def sync_pubmed(background_tasks: BackgroundTasks):
    job_id = f"sync_pubmed_{uuid4().hex[:8]}"
    data = {"id": job_id, "status": "queued", "source": "pubmed", "started_at": datetime.utcnow().isoformat()}
    supabase.table("sync_jobs").insert(data).execute()
    background_tasks.add_task(process_pubmed_data, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="PubMed sync started.")

@router.post("/sync/openfda", response_model=SyncJobResponse)
async def sync_openfda(background_tasks: BackgroundTasks):
    job_id = f"sync_openfda_{uuid4().hex[:8]}"
    data = {"id": job_id, "status": "queued", "source": "openfda", "started_at": datetime.utcnow().isoformat()}
    supabase.table("sync_jobs").insert(data).execute()
    background_tasks.add_task(openfda_service.sync_to_db, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="OpenFDA sync started.")

@router.post("/crawler/creative/run", response_model=SyncJobResponse)
async def run_creative_crawler(background_tasks: BackgroundTasks, search_term: Optional[str] = None):
    from app.services.creative_biolabs_crawler import creative_crawler
    job_id = f"crawl_creative_{uuid4().hex[:8]}"
    data = {"id": job_id, "status": "queued", "source": "creative_biolabs", "started_at": datetime.utcnow().isoformat()}
    supabase.table("sync_jobs").insert(data).execute()
    background_tasks.add_task(creative_crawler.run, search_term, 3, job_id)
    return SyncJobResponse(job_id=job_id, status="queued", message="Creative Biolabs crawler started.")

@router.get("/sync/{job_id}")
async def get_sync_status(job_id: str):
    job = await get_job_from_db(job_id)
    if not job: raise HTTPException(status_code=404, detail="Job not found")
    return job

@router.post("/sync/{job_id}/stop")
async def stop_sync_job(job_id: str):
    await update_job_status(job_id, cancel_requested=True)
    return {"message": "Stop request sent to DB."}

# --- Bulk Import & AI Refiner Endpoints ---

@router.post("/bulk/import", response_model=SyncJobResponse)
async def run_bulk_import(background_tasks: BackgroundTasks, max_studies: int = 5000):
    """ClinicalTrials.gov 전체 덤프에서 ADC 데이터 일괄 임포트"""
    from app.services.bulk_importer import BulkImporter
    
    job_id = f"bulk_import_{uuid4().hex[:8]}"
    data = {
        "id": job_id, 
        "status": "queued", 
        "source": "clinical_trials_bulk", 
        "started_at": datetime.utcnow().isoformat()
    }
    supabase.table("sync_jobs").insert(data).execute()
    
    importer = BulkImporter()
    background_tasks.add_task(importer.run_import, job_id, max_studies)
    
    return SyncJobResponse(
        job_id=job_id, 
        status="queued", 
        message=f"Bulk import started (max {max_studies} studies)."
    )

@router.post("/refiner/run", response_model=SyncJobResponse)
async def run_ai_refiner(background_tasks: BackgroundTasks, max_records: int = 50):
    """미정제 레코드 LLM 분석 및 SMILES 보강"""
    from app.services.ai_refiner import ai_refiner
    
    job_id = f"ai_refiner_{uuid4().hex[:8]}"
    data = {
        "id": job_id, 
        "status": "queued", 
        "source": "ai_refiner", 
        "started_at": datetime.utcnow().isoformat()
    }
    supabase.table("sync_jobs").insert(data).execute()
    
    background_tasks.add_task(ai_refiner.process_pending_records, job_id, max_records)
    
    return SyncJobResponse(
        job_id=job_id, 
        status="queued", 
        message=f"AI Refiner started (max {max_records} records)."
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
