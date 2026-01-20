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

from app.services.perplexity_service import enrich_drug_data
from app.core.config import settings
from app.core.supabase import supabase
from langchain_openai import ChatOpenAI
from langchain_core.prompts import ChatPromptTemplate
import json

router = APIRouter()

# Entrez 설정
Entrez.email = settings.ENTREZ_EMAIL

# [Cost Strategy] 데이터 처리용 가성비 모델
llm_refiner = ChatOpenAI(
    model=settings.FAST_LLM, # gpt-4o-mini
    temperature=0,
    api_key=settings.OPENAI_API_KEY
)

async def refine_drug_data_with_llm(title: str, description: str, raw_status: str, why_stopped: str) -> dict:
    """
    LLM을 사용하여 약물 이름 추출 및 임상 결과 분석 수행
    """
    system_prompt = """You are an expert Clinical Data Analyst. 
    Your task is to extract structured information from unstructured clinical trial text.
    
    1. **Drug Name Extraction:** 
       - Find the specific ADC code name (e.g., 'DS-8201', 'IMGN-632') or generic name.
       - If only 'Anti-HER2 ADC' is found, output 'Unknown'.
       - Do NOT invent names.
       
    2. **Outcome Analysis:**
       - Determine if the trial was a 'Success', 'Failure', or 'Ongoing'.
       - If 'Failure' (Terminated/Withdrawn), summarize the *scientific reason* (Toxicity, Lack of Efficacy) in 1 short sentence.
       
    Output JSON format:
    {
        "extracted_name": "DS-8201a",
        "outcome_type": "Failure", 
        "failure_reason": "Terminated due to Grade 3 interstitial lung disease."
    }
    """
    
    user_prompt = f"""
    Title: {title}
    Description: {description}
    Status: {raw_status}
    Why Stopped: {why_stopped}
    """
    
    prompt = ChatPromptTemplate.from_messages([
        ("system", system_prompt),
        ("user", user_prompt)
    ])
    
    try:
        chain = prompt | llm_refiner
        # JSON 모드 강제 (Structured Output)
        response = await chain.ainvoke({})
        
        content = response.content.strip()
        if "```json" in content:
            content = content.split("```json")[1].split("```")[0]
        elif "```" in content:
            content = content.split("```")[1].split("```")[0]
            
        return json.loads(content)
        
    except Exception as e:
        print(f"LLM Refine Error: {e}")
        return {"extracted_name": "Unknown", "outcome_type": "Unknown", "failure_reason": None}


# In-memory 동기화 상태 추적 (TODO: Redis로 대체)
sync_jobs: Dict[str, Dict[str, Any]] = {}


class SyncJobResponse(BaseModel):
    """동기화 작업 응답"""
    job_id: str
    status: str
    message: str


class SyncJobStatus(BaseModel):
    """동기화 작업 상태"""
    job_id: str
    status: str
    source: str
    records_found: int
    records_drafted: int
    started_at: str
    completed_at: Optional[str] = None
    errors: List[str] = []


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
async def fetch_clinical_trials(query: str = "ADC OR Antibody-Drug Conjugate"):
    """ClinicalTrials.gov API v2 호출"""
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.term": query,
        "pageSize": 20,
        "format": "json"
    }
    async with httpx.AsyncClient() as client:
        response = await client.get(url, params=params, timeout=30.0)
        response.raise_for_status()
        return response.json()


async def process_clinical_trials_data(job_id: str):
    """
    ClinicalTrials.gov 데이터 수집 및 처리
    """
    sync_jobs[job_id]["status"] = "running"
    
    try:
        # 1. Fetch Data
        data = await fetch_clinical_trials()
        studies = data.get("studies", [])
        
        sync_jobs[job_id]["records_found"] = len(studies)
        drafted = 0
        
        for study in studies:
            try:
                protocol = study.get("protocolSection", {})
                id_module = protocol.get("identificationModule", {})
                nct_id = id_module.get("nctId")
                
                # 중복 체크 (JSONB contains 사용)
                # existing = supabase.table("golden_set_library").select("id").eq("raw_data->>nct_id", nct_id).execute()
                existing = supabase.table("golden_set_library").select("id").contains("properties", {"nct_id": nct_id}).execute()
                
                if existing.data:
                    print(f"Skipping duplicate: {nct_id}")
                    continue

                # 기본 정보 추출
                conditions = protocol.get("conditionsModule", {}).get("conditions", [])
                interventions = protocol.get("armsInterventionsModule", {}).get("interventions", [])
                
                drug_name = "Unknown"
                for intervention in interventions:
                    if intervention.get("type") == "DRUG":
                        drug_name = intervention.get("name")
                        break
                
                record = {
                    "nct_id": nct_id,
                    "title": id_module.get("officialTitle"),
                    "conditions": conditions,
                    "drug_name": drug_name,
                    "phase": protocol.get("designModule", {}).get("phases", ["Unknown"])[0]
                }
                
                # 2. Perplexity Enrichment (Optional - 비용 절감을 위해 일부만 수행하거나 생략 가능)
                # enriched = await enrich_drug_data(record) 
                # 일단 Raw Data 위주로 저장
                
                # [LLM Injection] 데이터 정제 요청
                # 이름이 없거나, 상태가 'Terminated'인 경우 LLM에게 분석 요청
                should_refine = (
                    drug_name in ["Unknown", "Drug:"] or 
                    protocol.get("statusModule", {}).get("overallStatus") in ["TERMINATED", "WITHDRAWN", "SUSPENDED"]
                )
                
                refined_data = {}
                if should_refine:
                    refined_data = await refine_drug_data_with_llm(
                        title=id_module.get("officialTitle"),
                        description=str(protocol.get("descriptionModule", {})),
                        raw_status=protocol.get("statusModule", {}).get("overallStatus"),
                        why_stopped=protocol.get("statusModule", {}).get("whyStopped", "")
                    )
                
                # [Merge] LLM이 찾은 이름 적용
                final_name = drug_name
                if refined_data.get("extracted_name") and refined_data["extracted_name"] != "Unknown":
                    final_name = refined_data["extracted_name"]
                
                # 이름이 여전히 없으면 Fallback 규칙 적용
                if final_name == "Unknown":
                     final_name = f"Unnamed ADC ({nct_id})"

                # [Merge] LLM이 분석한 결과 적용
                outcome_type = refined_data.get("outcome_type", "Unknown")
                # LLM이 분석 안했으면 기본 매핑
                if outcome_type == "Unknown":
                    overall_status = protocol.get("statusModule", {}).get("overallStatus", "Unknown")
                    if overall_status in ["COMPLETED", "APPROVED"]:
                        outcome_type = "Success"
                    elif overall_status in ["TERMINATED", "WITHDRAWN", "SUSPENDED"]:
                        outcome_type = "Failure"
                    elif overall_status in ["RECRUITING", "ACTIVE_NOT_RECRUITING", "ENROLLING_BY_INVITATION"]:
                        outcome_type = "Ongoing"

                failure_reason = refined_data.get("failure_reason")
                if not failure_reason and outcome_type == "Failure":
                    failure_reason = protocol.get("statusModule", {}).get("whyStopped", "")

                new_entry = {
                    "name": final_name,
                    "category": "clinical_trial",
                    "description": f"[{outcome_type}] {record['title']}",
                    "properties": {
                        **record,
                        "ai_analysis": refined_data, # AI 분석 원본 보관
                        "overall_status": protocol.get("statusModule", {}).get("overallStatus"),
                        "why_stopped": protocol.get("statusModule", {}).get("whyStopped")
                    },
                    "status": "draft",
                    "enrichment_source": "clinical_trials_ai_refined" if should_refine else "clinical_trials",
                    "raw_data": study, # 전체 데이터 저장
                    "outcome_type": outcome_type,
                    "failure_reason": failure_reason
                }

                supabase.table("golden_set_library").insert(new_entry).execute()
                
                drafted += 1
                sync_jobs[job_id]["records_drafted"] = drafted
                
                # Rate Limiting
                await asyncio.sleep(0.5)
                
            except Exception as e:
                print(f"Error processing {nct_id}: {e}")
                sync_jobs[job_id]["errors"].append(f"{nct_id}: {str(e)}")
        
        sync_jobs[job_id]["status"] = "completed"
        sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        
    except Exception as e:
        sync_jobs[job_id]["status"] = "failed"
        sync_jobs[job_id]["errors"].append(str(e))
        
        # Log to DB
        try:
            supabase.table("data_sync_logs").insert({
                "source_id": "clinical_trials",
                "status": "failed",
                "records_synced": 0,
                "records_drafted": 0,
                "error_message": str(e)[:1000]
            }).execute()
        except:
            pass


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=1, min=2, max=10))
def fetch_pubmed_ids(term: str = "Antibody-Drug Conjugate", max_results: int = 20):
    """PubMed ID 검색"""
    handle = Entrez.esearch(db="pubmed", term=term, retmax=max_results)
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
    """
    PubMed 데이터 수집
    """
    sync_jobs[job_id]["status"] = "running"
    
    try:
        # 1. Search IDs
        id_list = fetch_pubmed_ids()
        sync_jobs[job_id]["records_found"] = len(id_list)
        
        if not id_list:
            sync_jobs[job_id]["status"] = "completed"
            return

        # 2. Fetch Details
        # Entrez는 동기 함수이므로, 별도 스레드나 비동기 래퍼 없이 호출하면 블로킹됨.
        # 여기서는 간단히 호출 (BackgroundTasks 내에서 실행되므로 메인 스레드 차단 안함)
        papers = fetch_pubmed_details(id_list)
        
        drafted = 0
        if 'PubmedArticle' in papers:
            for article in papers['PubmedArticle']:
                try:
                    medline = article['MedlineCitation']
                    article_data = medline['Article']
                    
                    pmid = str(medline['PMID'])
                    title = article_data.get('ArticleTitle', 'No Title')
                    abstract_list = article_data.get('Abstract', {}).get('AbstractText', [])
                    abstract = " ".join(abstract_list) if abstract_list else "No Abstract"
                    
                    # 중복 체크 (Title은 text 컬럼이 없으므로 properties나 raw_data 확인 필요하지만, 
                    # knowledge_base는 title 컬럼이 있음! -> eq 사용 가능)
                    # 단, source_type과 title 모두 일치해야 함
                    existing = supabase.table("knowledge_base").select("id").eq("source_type", "PubMed").eq("title", title).execute()
                    if existing.data:
                        print(f"Skipping duplicate PubMed: {pmid}")
                        continue
                        
                    # Knowledge Base에 저장
                    new_kb = {
                        "source_type": "PubMed",
                        "title": title,
                        "summary": abstract[:500] + "...", # 요약은 나중에 AI로
                        "content": abstract,
                        "relevance_score": 0.0, # 초기값
                        "source_tier": 1, # PubMed는 Tier 1
                        "rag_status": "pending"
                    }
                    
                    supabase.table("knowledge_base").insert(new_kb).execute()
                    
                    drafted += 1
                    sync_jobs[job_id]["records_drafted"] = drafted
                    
                except Exception as e:
                    sync_jobs[job_id]["errors"].append(f"PMID {pmid}: {str(e)}")

        sync_jobs[job_id]["status"] = "completed"
        sync_jobs[job_id]["completed_at"] = datetime.utcnow().isoformat()
        
    except Exception as e:
        sync_jobs[job_id]["status"] = "failed"
        sync_jobs[job_id]["errors"].append(str(e))


@router.post("/sync/clinical", response_model=SyncJobResponse)
async def sync_clinical_trials(background_tasks: BackgroundTasks):
    """ClinicalTrials.gov 데이터 동기화"""
    job_id = f"sync_clinical_{uuid4().hex[:8]}"
    
    sync_jobs[job_id] = {
        "status": "queued",
        "source": "clinical_trials",
        "records_found": 0,
        "records_drafted": 0,
        "started_at": datetime.utcnow().isoformat(),
        "completed_at": None,
        "errors": []
    }
    
    background_tasks.add_task(process_clinical_trials_data, job_id)
    
    return SyncJobResponse(
        job_id=job_id,
        status="queued",
        message="ClinicalTrials.gov sync started."
    )


@router.post("/sync/pubmed", response_model=SyncJobResponse)
async def sync_pubmed(background_tasks: BackgroundTasks):
    """PubMed 데이터 동기화"""
    job_id = f"sync_pubmed_{uuid4().hex[:8]}"
    
    sync_jobs[job_id] = {
        "status": "queued",
        "source": "pubmed",
        "records_found": 0,
        "records_drafted": 0,
        "started_at": datetime.utcnow().isoformat(),
        "completed_at": None,
        "errors": []
    }
    
    background_tasks.add_task(process_pubmed_data, job_id)
    
    return SyncJobResponse(
        job_id=job_id,
        status="queued",
        message="PubMed sync started."
    )


@router.get("/sync/{job_id}", response_model=SyncJobStatus)
async def get_sync_status(job_id: str):
    """동기화 작업 상태 조회"""
    if job_id not in sync_jobs:
        raise HTTPException(status_code=404, detail="Sync job not found")
    
    job = sync_jobs[job_id]
    return SyncJobStatus(job_id=job_id, **job)


@router.get("/health")
async def scheduler_health():
    """스케줄러 API 헬스 체크"""
    return {
        "status": "healthy",
        "active_jobs": len([j for j in sync_jobs.values() if j["status"] == "running"])
    }
