"""
PubMed Knowledge Base Service
golden_set_library 기반 PubMed 논문 자동 수집 및 Gemini AI 분석
내일 오전 9시까지 약 2,000개 약물에 대한 논문 수집 목표
"""
import asyncio
import logging
import json
from typing import List, Dict, Any, Optional
from datetime import datetime, timedelta
from Bio import Entrez

import google.generativeai as genai
from json_repair import repair_json

from app.core.config import settings
from app.core.supabase import supabase
from app.services.cost_tracker import cost_tracker
from app.services.chemical_mapper import chemical_mapper
from app.services.rag_service import rag_service

logger = logging.getLogger(__name__)

# Entrez 설정
Entrez.email = settings.NCBI_EMAIL
if settings.NCBI_API_KEY:
    Entrez.api_key = settings.NCBI_API_KEY


class PubMedKnowledgeService:
    """golden_set_library 기반 PubMed 논문 수집 및 AI 분석"""
    
    # Rate Limiting 설정
    # API Key 있으면 10 req/sec, 없으면 3 req/sec
    REQUESTS_PER_SECOND = 10 if settings.NCBI_API_KEY else 3
    REQUEST_DELAY = 1.0 / REQUESTS_PER_SECOND
    MAX_PAPERS_PER_DRUG = 10  # 골든셋 전수 수집을 위해 상향 조정
    YEARS_BACK = 10 # 최근 10년으로 확장
    
    # Gemini Safety Settings (의학 용어 차단 해제)
    SAFETY_SETTINGS = [
        {"category": "HARM_CATEGORY_HATE_SPEECH", "threshold": "BLOCK_NONE"},
        {"category": "HARM_CATEGORY_HARASSMENT", "threshold": "BLOCK_NONE"},
        {"category": "HARM_CATEGORY_SEXUALLY_EXPLICIT", "threshold": "BLOCK_NONE"},
        {"category": "HARM_CATEGORY_DANGEROUS_CONTENT", "threshold": "BLOCK_NONE"}
    ]
    
    def __init__(self):
        self.processed_count = 0
        self.error_count = 0
        self.skipped_count = 0
        self.failed_drugs = []  # 실패 로그
    
    async def get_target_drugs(self, limit: Optional[int] = None) -> List[Dict]:
        """
        OpenFDA 외 enrichment_source를 가진 약물 리스트 추출
        또는 데이터가 부족한 약물 추출
        """
        try:
            query = supabase.table("golden_set_library")\
                .select("id, name, properties")\
                .neq("enrichment_source", "open_fda_api")
            
            if limit:
                query = query.limit(limit)
            
            response = query.execute()
            
            drugs = []
            for record in response.data:
                name = record.get("name")
                if name and name != "Unknown":
                    props = record.get("properties", {}) or {}
                    drugs.append({
                        "id": record["id"],
                        "name": name,
                        "generic_name": props.get("generic_name")
                    })
            
            logger.info(f"📋 Found {len(drugs)} target drugs from golden_set_library")
            return drugs
            
        except Exception as e:
            logger.error(f"Failed to fetch target drugs: {e}")
            return []
    
    def build_search_query(self, drug_name: str, generic_name: str = None) -> str:
        """
        PubMed 검색 쿼리 생성 - MeSH Term 활용 고도화
        형식: "{Drug Name}" AND (MeSH Terms OR "ADC")
        """
        # 약물명 정리 (특수문자 제거)
        clean_name = drug_name.replace('"', '').strip()
        
        # 기본 쿼리 (제목/초록)
        base_query = f'("{clean_name}"[Title/Abstract])'
        
        # MeSH Terms 및 동의어 활용 (수집 범위 확장 - Mega-Net Strategy)
        # ADC, 항체, 이중항체, 접합체, 그리고 주요 접미사 패턴까지 모두 포함
        adc_context = (
            '('
            # 1. MeSH & General Terms
            '"Antibody-Drug Conjugates"[MeSH] OR "Antibody-Drug Conjugate" OR "ADC" OR '
            '"Immunoconjugates"[MeSH] OR "Immunoconjugate" OR '
            '"Antibodies, Monoclonal"[MeSH] OR "Monoclonal Antibody" OR '
            '"Bispecific Antibodies"[MeSH] OR "Bispecific Antibody" OR "bsAb" OR '
            '"Antibodies, Monoclonal, Conjugated"[MeSH] OR '
            '"Checkpoint Inhibitor" OR '
            
            # 2. Suffix Patterns (Wildcards)
            '*-mab OR *-vedotin OR *-deruxtecan OR *-emtansine OR '
            '*-ozogamicin OR *-mafodotin OR *-tesirine OR *-voralatide OR *-stansine'
            ')'
        )
        
        # 최종 쿼리
        query = f'{base_query} AND {adc_context}'
        
        return query
    
    async def search_general_discovery(self, max_results: int = 10000) -> List[Dict]:
        """
        [전수 수집용] Yearly Chunking Strategy + Mega-Net Query
        2010년부터 2026년까지 연도별로 분할하여 대량 수집 (목표: 1만 건)
        """
        loop = asyncio.get_event_loop()
        
        # 1. Mega-Net Query Expansion
        # *-mab, *-tin, *-can 등 포괄적 접미사 및 ADC 관련 키워드 대폭 확장
        query_terms = [
            # Core ADC Terms
            '"Antibody-Drug Conjugates"[MeSH]', '"Antibody-Drug Conjugate"', '"ADC"',
            '"Immunoconjugates"[MeSH]', '"Immunoconjugate"',
            '"Antibodies, Monoclonal, Conjugated"[MeSH]',
            
            # Antibody Terms
            '"Antibodies, Monoclonal"[MeSH]', '"Monoclonal Antibody"',
            '"Bispecific Antibodies"[MeSH]', '"Bispecific Antibody"', '"bsAb"',
            '"Checkpoint Inhibitor"',
            
            # Suffix Wildcards (Text Word) - 검색 범위 5배 확장
            '*-mab', '*-tin', '*-can', '*-tecan', '*-tansine', 
            '*-vedotin', '*-deruxtecan', '*-emtansine', 
            '*-ozogamicin', '*-mafodotin', '*-tesirine', 
            '*-voralatide', '*-stansine'
        ]
        
        # Combine with OR
        main_query = " OR ".join(query_terms)
        
        # Filter (Clinical Focus but inclusive enough)
        # 단순히 Clinical Trial만 찾지 않고, 새로운 전임상/연구 결과도 포함하도록 필터 완화
        # 다만 너무 noise가 많으면 조정 필요. 현재는 Discovery 모드이므로 넓게 잡음.
        full_query = f'({main_query})' 
        
        logger.info(f"📡 Mega-Net Query Strategy Activated")
        logger.info(f"   Query: {full_query[:200]}...")

        all_articles = []
        current_year = datetime.now().year
        start_year = 2010
        
        # 2. Yearly Chunking Loop
        try:
            for year in range(current_year, start_year - 1, -1):
                if len(all_articles) >= max_results:
                    break
                
                logger.info(f"   📅 Fetching Year: {year}...")
                
                def do_search_year(target_year):
                    handle = Entrez.esearch(
                        db="pubmed",
                        term=full_query,
                        retmax=5000, # 연도별 최대치 (넉넉하게)
                        sort="relevance",
                        datetype="pdat",
                        mindate=f"{target_year}/01/01",
                        maxdate=f"{target_year}/12/31"
                    )
                    record = Entrez.read(handle)
                    handle.close()
                    return record.get("IdList", [])
                
                id_list = await loop.run_in_executor(None, lambda: do_search_year(year))
                logger.info(f"      -> Found {len(id_list)} IDs in {year}")
                
                if not id_list:
                    continue
                
                # Fetch details in sub-batches (Entrez limit handling)
                chunk_size = 50 # 200 -> 50으로 축소 (반응성 향상)
                for i in range(0, len(id_list), chunk_size):
                    if len(all_articles) >= max_results:
                        break
                        
                    chunk_ids = id_list[i : i + chunk_size]
                    
                    # 통신 시작 전 로그
                    logger.info(f"      ⏳ Fetching details for {len(chunk_ids)} papers (Batch {i // chunk_size + 1})...")
                    
                    articles = await self._fetch_details(chunk_ids)
                    
                    # Add processing context
                    for a in articles:
                        a["drug_name"] = "Discovery" # Placeholder
                    
                    all_articles.extend(articles)
                    logger.info(f"      ✅ Fetched {len(articles)} details (Total: {len(all_articles)})")
                    
                    await asyncio.sleep(self.REQUEST_DELAY) # Rate limit compliance
                
            return all_articles[:max_results]
            
        except Exception as e:
            logger.error(f"General discovery search error: {e}")
            return all_articles # Return what we have so far

    async def _fetch_details(self, id_list: List[str]) -> List[Dict]:
        """PMID 리스트로 상세 정보 조회 (공통 함수)"""
        loop = asyncio.get_event_loop()
        if not id_list: return []
        
        try:
            await asyncio.sleep(self.REQUEST_DELAY)
            def fetch_details():
                handle = Entrez.efetch(
                    db="pubmed",
                    id=",".join(id_list),
                    retmode="xml"
                )
                records = Entrez.read(handle)
                handle.close()
                return records
            
            records = await loop.run_in_executor(None, fetch_details)
            
            articles = []
            for article in records.get("PubmedArticle", []):
                try:
                    medline = article["MedlineCitation"]
                    article_data = medline["Article"]
                    
                    title = str(article_data.get("ArticleTitle", "No Title"))
                    abstract_texts = article_data.get("Abstract", {}).get("AbstractText", [])
                    if isinstance(abstract_texts, list):
                        abstract = " ".join([str(t) for t in abstract_texts])
                    else:
                        abstract = str(abstract_texts) if abstract_texts else ""
                    
                    pmid = str(medline.get("PMID", ""))
                    journal = article_data.get("Journal", {}).get("Title", "")
                    pub_date = article_data.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
                    year = pub_date.get("Year", "")
                    
                    articles.append({
                        "pmid": pmid,
                        "title": title,
                        "abstract": abstract,
                        "journal": journal,
                        "year": year,
                        "drug_name": "Discovery" # Placeholder
                    })
                except Exception as e:
                    logger.warning(f"Failed to parse article: {e}")
            return articles
        except Exception as e:
            logger.error(f"Fetch details error: {e}")
            return []

    async def search_pubmed_for_drug(
        self, 
        drug_name: str, 
        generic_name: str = None,
        max_results: int = 5
    ) -> List[Dict]:
        """약물별 PubMed 논문 검색"""
        loop = asyncio.get_event_loop()
        
        # 1차 검색
        query = self.build_search_query(drug_name)
        
        try:
            end_date = datetime.now()
            start_date = end_date - timedelta(days=365 * self.YEARS_BACK)
            
            def do_search():
                handle = Entrez.esearch(
                    db="pubmed",
                    term=query,
                    retmax=max_results,
                    sort="relevance",
                    datetype="pdat",
                    mindate=start_date.strftime("%Y/%m/%d"),
                    maxdate=end_date.strftime("%Y/%m/%d")
                )
                record = Entrez.read(handle)
                handle.close()
                return record.get("IdList", [])
            
            id_list = await loop.run_in_executor(None, do_search)
            
            # Fallback
            if not id_list and generic_name and generic_name != drug_name:
                logger.info(f"🔄 Fallback: Searching with generic name '{generic_name}'")
                query = self.build_search_query(generic_name)
                id_list = await loop.run_in_executor(None, do_search)
            
            if not id_list:
                return []
            
            articles = await self._fetch_details(id_list)
            for a in articles: a["drug_name"] = drug_name # Restore specific drug name context
            return articles
            
        except Exception as e:
            logger.error(f"PubMed search error for '{drug_name}': {e}")
            return []
    
    async def analyze_with_gemini(self, abstract: str, title: str, drug_name: str) -> Dict[str, Any]:
        """
        Gemini 2.0 Flash 분석 - 고도화 (Structured Output for ADC Platform)
        4단 분리 추출 (Target / Antibody / Linker / Payload) 및 임상 수치 정형화 강화
        """
        # 1. Regex Extraction (Safety Net)
        import re
        
        # Extract NCT IDs (NCT + 8 digits)
        nct_ids = re.findall(r"NCT\d{8}", abstract + " " + title)
        nct_id_val = nct_ids[0] if nct_ids else None
        
        system_prompt = """You are an expert ADC (Antibody-Drug Conjugate) researcher.
Analyze the abstract and extract structured data into a specific JSON format.

CRITICAL INSTRUCTIONS:
1. **4-Part Decomposition**: You MUST attempt to separate the ADC into:
   - Target Antigen (e.g., HER2, TROP2, CD19)
   - Antibody (e.g., Trastuzumab, Sacituzumab)
   - Linker (e.g., GGFG peptide, hydrazone, mc-vc)
   - Payload (e.g., DXd, SN-38, MMAE)
   *If a specific part is not mentioned, use null.*

2. **Clinical Metrics**: Extract explicit numbers for:
   - ORR (Objective Response Rate) in %
   - PFS (Progression-Free Survival) in months
   - OS (Overall Survival) in months
   - DoR (Duration of Response) in months

3. **Drug Name**: Identify the MAIN ADC/Drug code name (e.g., DS-8201, ABBV-400).

Output MUST be a JSON object with this EXACT structure:
{
  "drug_structure": {
    "drug_name": "Main drug name",
    "target": "Target Antigen or null",
    "antibody": "Antibody part or null",
    "linker": "Linker part or null",
    "payload": "Payload part or null"
  },
  "clinical_metrics": {
    "nct_id": "Clinical Trial ID (e.g., NCT03529110) or null",
    "orr_pct": Float or null,
    "pfs_months": Float or null,
    "os_months": Float or null,
    "dor_months": Float or null
  },
  "analysis": {
    "outcome_type": "Success", "Failure", "Ongoing", or "Unknown",
    "reasoning": "Brief explanation of the trial result or study conclusion",
    "is_golden_set_candidate": boolean (true if it contains significant clinical data),
    "ai_confidence_score": Float (0.0-1.0),
    "summary": "Concise 3-sentence summary highlighting efficacy and safety"
  }
}

IMPORTANT: Return ONLY raw JSON. No markdown."""

        full_prompt = f"""{system_prompt}

Drug of Interest Context: {drug_name}
Title: {title}
Abstract: {abstract[:3500]}"""

        try:
            genai.configure(api_key=settings.GOOGLE_API_KEY)
            model_id = settings.GEMINI_MODEL_ID or 'gemini-2.0-flash'
            model = genai.GenerativeModel(model_id)
            
            loop = asyncio.get_event_loop()
            response = await loop.run_in_executor(
                None, 
                lambda: model.generate_content(
                    full_prompt,
                    safety_settings=self.SAFETY_SETTINGS
                )
            )
            
            content = response.text.strip()
            usage = response.usage_metadata
            await cost_tracker.track_usage(model_id, usage.prompt_token_count, usage.candidates_token_count)
            
            try:
                repaired = repair_json(content)
                result = json.loads(repaired)
            except:
                if "```json" in content:
                    content = content.split("```json")[1].split("```")[0]
                elif "```" in content:
                    content = content.split("```")[1].split("```")[0]
                result = json.loads(content.strip())
            
            # --- Post-Processing & Validation ---
            
            # 1. Force Inject Regex NCT ID (Priority: Regex > AI)
            if nct_id_val:
                if "clinical_metrics" not in result: result["clinical_metrics"] = {}
                result["clinical_metrics"]["nct_id"] = nct_id_val
                
            # 2. Chemical Mapper Integration
            extracted_drug_name = result.get("drug_structure", {}).get("drug_name")
            payload_name = result.get("drug_structure", {}).get("payload")
            
            # Only run mapper if we actually found a payload name
            if payload_name and payload_name.lower() not in ["null", "none"]:
                try:
                    mapped_data = await chemical_mapper.enrich_with_commercial_data(payload_name)
                    if mapped_data.get("payload_smiles"):
                        result["drug_structure"]["payload_smiles"] = mapped_data["payload_smiles"]
                except Exception as map_err:
                    logger.warning(f"Mapper failed for payload '{payload_name}': {map_err}")

            return result
                
        except Exception as e:
            logger.error(f"Gemini analysis error: {e}")
            return {
                "analysis": {
                    "summary": "Analysis failed",
                    "ai_confidence_score": 0.0,
                    "reasoning": f"Error: {str(e)}"
                }
            }
    
    async def is_pmid_duplicate(self, pmid: str) -> Optional[str]:
        """PMID 기준 중복 검사 - ID 반환"""
        if not pmid: return None
        try:
            # 2. Summary text search for "PMID: {pmid}"
            # Note: A separate 'pmid' column would be better, but sticking to current schema
            existing = supabase.table("knowledge_base")\
                .select("id")\
                .ilike("summary", f"%PMID: {pmid}%")\
                .limit(1)\
                .execute()
            
            if existing.data:
                return existing.data[0]["id"]
            return None
        except:
            return None
    
    async def save_to_knowledge_base(
        self, 
        article: Dict, 
        analysis: Dict,
        drug_id: Optional[str] = None
    ) -> bool:
        """
        1. PMID Check (Dedup & Upsert)
        2. Generate Embedding
        3. Save (Insert or Update) to knowledge_base
        4. Update golden_set_library (if linked)
        """
        try:
            # Extract fields from nested structure
            an = analysis.get("analysis", {})
            ds = analysis.get("drug_structure", {})
            cm = analysis.get("clinical_metrics", {})
            
            pmid = article.get("pmid")
            summary_text = f"PMID: {pmid} | " + an.get("summary", "No summary provided.")
            reasoning = an.get("reasoning", "No reasoning provided.")
            
            # 1. Dedup / Upsert Check
            existing_id = await self.is_pmid_duplicate(pmid)
            
            if existing_id:
                logger.info(f"🔄 Duplicate PMID found ({pmid}). Updating properties (Upsert)...")
                
                # Upsert Logic: Only update analysis/properties
                update_data = {
                    "properties": analysis,
                    "ai_reasoning": reasoning[:2000],
                    "relevance_score": an.get("ai_confidence_score", 0.0),
                    # Optionally update summary if needed, but keeping it stable might be better
                }
                supabase.table("knowledge_base").update(update_data).eq("id", existing_id).execute()
                return True

            # 3. Generate Embedding (Gemini 768) - Only for new records
            embedding_text = f"{article['title']} {summary_text} {reasoning}"
            embedding = await rag_service.generate_embedding(embedding_text)

            new_kb = {
                "source_type": "PubMed",
                "title": article["title"][:500],
                "content": article.get("abstract", "No abstract available."),
                "summary": summary_text[:1000],
                "relevance_score": an.get("ai_confidence_score", 0.0),
                "source_tier": 1,
                "ai_reasoning": reasoning[:2000],
                "rag_status": "completed",
                "embedding": embedding,
                "properties": analysis # Save the FULL nested JSON structure here
            }
            
            supabase.table("knowledge_base").insert(new_kb).execute()
            
            # 4. Update Golden Set Library (Linkage)
            if drug_id and an.get("ai_confidence_score", 0) >= 0.7:
                await self._update_golden_set(drug_id, analysis)
            
            return True
            
        except Exception as e:
            logger.error(f"Failed to save to knowledge_base: {e}")
            return False

    async def _update_golden_set(self, drug_id: str, analysis: Dict):
        """골든셋 라이브러리 업데이트 (별도 함수 분리)"""
        try:
            curr_res = supabase.table("golden_set_library").select("properties, outcome_type, failure_reason").eq("id", drug_id).single().execute()
            current_props = curr_res.data.get("properties", {}) or {}
            
            # Extract flat fields for backward compatibility in golden_set_library, or store nested?
            # Storing nested is cleaner, but let's flatten key metrics to top level of properties for easy UI access if needed
            ds = analysis.get("drug_structure", {})
            cm = analysis.get("clinical_metrics", {})
            an = analysis.get("analysis", {})
            
            updates = {}
            
            # Update specific fields if they exist
            if cm.get("orr_pct"): current_props["orr"] = f"{cm['orr_pct']}%"
            if cm.get("pfs_months"): current_props["pfs"] = f"{cm['pfs_months']} mo"
            if cm.get("os_months"): current_props["os"] = f"{cm['os_months']} mo"
            
            if ds.get("linker"): current_props["linker_type"] = ds["linker"]
            if ds.get("payload"): current_props["payload_name"] = ds["payload"]
            if ds.get("target"): current_props["target_antigen"] = ds["target"]
            
            if an.get("outcome_type"): updates["outcome_type"] = an["outcome_type"]
            if an.get("reasoning"): updates["failure_reason"] = an["reasoning"]
            
            updates["properties"] = current_props
            updates["updated_at"] = datetime.utcnow().isoformat()
            
            supabase.table("golden_set_library").update(updates).eq("id", drug_id).execute()
            logger.info(f"✨ Updated Golden Set properties for drug {drug_id}")
            
        except Exception as e:
            logger.error(f"Failed to update golden_set_library: {e}")

    async def process_single_drug(self, drug: Dict, job_id: Optional[str] = None) -> int:
        """단일 약물에 대한 PubMed 논문 수집 및 분석"""
        drug_name = drug["name"]
        generic_name = drug.get("generic_name")
        drug_id = drug.get("id")
        
        saved_count = 0
        
        try:
            # 1. PubMed 검색
            articles = await self.search_pubmed_for_drug(
                drug_name, 
                generic_name, 
                max_results=self.MAX_PAPERS_PER_DRUG
            )
            
            if not articles:
                logger.debug(f"⏩ No articles found for: {drug_name}")
                return 0
            
            # 2. 각 논문 처리
            for article in articles:
                # 중복 체크 (PMID) - 상위 레벨에서 한 번 더 체크해도 좋지만 save_to_knowledge_base 내부에서 처리
                
                # 초록이 없으면 스킵
                if not article.get("abstract"):
                    continue
                
                # 3. AI 분석
                analysis = await self.analyze_with_gemini(
                    article["abstract"],
                    article["title"],
                    drug_name
                )
                
                # 4. 저장
                if await self.save_to_knowledge_base(article, analysis, drug_id):
                    saved_count += 1
                    self.processed_count += 1
                    logger.info(f"✅ Saved: {article['title'][:50]}... (Score: {analysis.get('relevance_score', 0):.2f})")
                
                # Rate limiting between AI calls
                await asyncio.sleep(0.3)
            
            return saved_count
            
        except Exception as e:
            self.error_count += 1
            self.failed_drugs.append({"drug_name": drug_name, "error": str(e)})
            logger.error(f"❌ Error processing drug '{drug_name}': {e}")
            return 0
    
    async def run_batch(
        self, 
        job_id: Optional[str] = None, 
        batch_size: int = 100,
        mode: str = "incremental"
    ) -> Dict[str, Any]:
        """
        배치 실행
        mode: 'incremental' (신규만), 'full' (전체), 'discovery' (광범위 검색)
        """
        from app.api.scheduler import update_job_status, is_cancelled
        
        if job_id:
            await update_job_status(job_id, status="running")
        
        logger.info(f"🚀 [PubMed Knowledge] Starting batch (size: {batch_size}, mode: {mode})")
        
        try:
            # Mode: Discovery (일반 검색)
            if mode == "discovery":
                logger.info(f"📥 Fetching details for found IDs...")
                articles = await self.search_general_discovery(max_results=batch_size)
                
                if not articles:
                    logger.info("ℹ️ No articles found to process.")
                    return {"status": "completed", "processed": 0}
                
                total_articles = len(articles)
                logger.info(f"⚡ Starting AI Analysis for {total_articles} articles...")
                
                for idx, article in enumerate(articles):
                    if job_id and await is_cancelled(job_id): 
                        logger.warning("🛑 Job cancelled by user.")
                        break
                    
                    # 진행 상황 로그 (CLI용)
                    if (idx + 1) % 10 == 0 or idx == 0:
                        logger.info(f"🔍 Analyzing {idx + 1}/{total_articles} ({(idx + 1)/total_articles*100:.1f}%)...")

                    analysis = await self.analyze_with_gemini(article["abstract"], article["title"], "General Discovery")
                    
                    if await self.save_to_knowledge_base(article, analysis):
                        self.processed_count += 1
                        logger.info(f"   ✅ Saved: {article['title'][:60]}... (Score: {analysis.get('analysis', {}).get('ai_confidence_score', 0):.2f})")
                    
                    await asyncio.sleep(0.3)
                    
                    if job_id and (idx + 1) % 5 == 0:
                        await update_job_status(job_id, records_drafted=self.processed_count)
                
                return {"status": "completed", "papers_saved": self.processed_count}

            # Mode: Standard (Drug list based)
            limit = batch_size if mode == "incremental" else None
            drugs = await self.get_target_drugs(limit)
            
            if not drugs:
                if job_id: await update_job_status(job_id, status="completed", message="No drugs")
                return {"status": "completed", "processed": 0}
            
            if job_id: await update_job_status(job_id, records_found=len(drugs))
            
            total_saved = 0
            for idx, drug in enumerate(drugs):
                if job_id and await is_cancelled(job_id): break
                
                saved = await self.process_single_drug(drug, job_id)
                total_saved += saved
                
                if job_id and (idx + 1) % 10 == 0:
                    await update_job_status(
                        job_id, 
                        records_drafted=self.processed_count,
                        message=f"Processing {idx + 1}/{len(drugs)} drugs..."
                    )
                await asyncio.sleep(self.REQUEST_DELAY)
            
            result = {
                "status": "completed",
                "total_drugs": len(drugs),
                "papers_saved": self.processed_count,
                "errors": self.error_count
            }
            
            if job_id:
                await update_job_status(
                    job_id,
                    status="completed",
                    records_drafted=self.processed_count,
                    message=f"Saved {self.processed_count} papers"
                )
            
            return result
            
        except Exception as e:
            logger.error(f"Batch processing error: {e}")
            if job_id: await update_job_status(job_id, status="failed", errors=[str(e)])
            return {"status": "failed", "error": str(e)}


# 싱글톤 인스턴스
pubmed_knowledge_service = PubMedKnowledgeService()