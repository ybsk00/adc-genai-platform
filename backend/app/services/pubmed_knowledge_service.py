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
    MAX_PAPERS_PER_DRUG = 5
    YEARS_BACK = 5  # 최근 5년 논문만
    
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
        PubMed 검색 쿼리 생성
        형식: "{Drug Name}" AND ("antibody-drug conjugate" OR "ADC")
        """
        # 약물명 정리 (특수문자 제거)
        clean_name = drug_name.replace('"', '').strip()
        
        # 기본 쿼리
        base_query = f'("{clean_name}"[Title/Abstract])'
        
        # ADC 컨텍스트 추가
        adc_context = '("antibody-drug conjugate"[Title/Abstract] OR "ADC"[Title/Abstract] OR "immunoconjugate"[Title/Abstract])'
        
        # 최종 쿼리
        query = f'{base_query} AND {adc_context}'
        
        return query
    
    async def search_pubmed_for_drug(
        self, 
        drug_name: str, 
        generic_name: str = None,
        max_results: int = 5
    ) -> List[Dict]:
        """
        약물별 PubMed 논문 검색
        필터: 최근 5년, English, Has Abstract
        Fallback: 브랜드명 실패 시 generic_name으로 재검색
        """
        loop = asyncio.get_event_loop()
        
        # 1차 검색: 브랜드명
        query = self.build_search_query(drug_name)
        
        try:
            # 날짜 범위 설정 (최근 5년)
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
            
            # Fallback: 결과 없으면 generic_name으로 재검색
            if not id_list and generic_name and generic_name != drug_name:
                logger.info(f"🔄 Fallback: Searching with generic name '{generic_name}'")
                query = self.build_search_query(generic_name)
                id_list = await loop.run_in_executor(None, do_search)
            
            if not id_list:
                return []
            
            # Rate limiting
            await asyncio.sleep(self.REQUEST_DELAY)
            
            # 상세 정보 조회
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
                    
                    # 제목 추출
                    title = str(article_data.get("ArticleTitle", "No Title"))
                    
                    # 초록 추출
                    abstract_texts = article_data.get("Abstract", {}).get("AbstractText", [])
                    if isinstance(abstract_texts, list):
                        abstract = " ".join([str(t) for t in abstract_texts])
                    else:
                        abstract = str(abstract_texts) if abstract_texts else ""
                    
                    # PMID 추출
                    pmid = str(medline.get("PMID", ""))
                    
                    # 저널 정보
                    journal = article_data.get("Journal", {}).get("Title", "")
                    
                    # 발행일
                    pub_date = article_data.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
                    year = pub_date.get("Year", "")
                    
                    articles.append({
                        "pmid": pmid,
                        "title": title,
                        "abstract": abstract,
                        "journal": journal,
                        "year": year,
                        "drug_name": drug_name
                    })
                except Exception as e:
                    logger.warning(f"Failed to parse article: {e}")
            
            return articles
            
        except Exception as e:
            logger.error(f"PubMed search error for '{drug_name}': {e}")
            return []
    
    async def analyze_with_gemini(self, abstract: str, title: str, drug_name: str) -> Dict[str, Any]:
        """
        Gemini 1.5 Flash로 논문 분석
        - Target 추출 (HER2, TROP2, CD19 등)
        - Indication 추출 (암종)
        - Summary (3문장 이내)
        - Relevance Score (0.0~1.0)
        - AI Reasoning (판단 근거)
        """
        system_prompt = """You are an expert ADC (Antibody-Drug Conjugate) researcher analyzing scientific publications.

Analyze the provided scientific abstract and extract structured information.

Output MUST be a JSON object with these exact fields:
1. "target": Molecular target(s) mentioned (e.g., "HER2", "TROP2", "CD19", "Nectin-4"). Return null if not found.
2. "indication": Cancer type or disease indication (e.g., "breast cancer", "acute lymphoblastic leukemia"). Return null if not found.
3. "summary": A concise 3-sentence summary focusing on:
   - Clinical stage (Phase 1, 2, 3 if mentioned)
   - Key efficacy results (ORR, PFS, OS if available)
   - Main conclusions
4. "relevance_score": Float between 0.0 and 1.0 for ADC design relevance:
   - 1.0: Directly about ADC clinical trials, efficacy, or safety
   - 0.7-0.9: ADC mechanism, linker, or payload studies
   - 0.4-0.6: General antibody or conjugate research
   - 0.1-0.3: Tangentially related
5. "ai_reasoning": One sentence explaining why this paper is important for ADC research.

IMPORTANT: Return ONLY raw JSON. Do not use markdown formatting."""

        full_prompt = f"""{system_prompt}

Drug of Interest: {drug_name}
Title: {title}
Abstract: {abstract[:3000]}"""  # 토큰 제한을 위해 초록 3000자로 제한

        try:
            genai.configure(api_key=settings.GOOGLE_API_KEY)
            model = genai.GenerativeModel('gemini-1.5-flash')  # 지시서대로 1.5 Flash 사용
            
            loop = asyncio.get_event_loop()
            response = await loop.run_in_executor(
                None, 
                lambda: model.generate_content(full_prompt)
            )
            
            content = response.text.strip()
            
            # 비용 추적
            usage = response.usage_metadata
            await cost_tracker.track_usage(
                "gemini-1.5-flash",
                usage.prompt_token_count,
                usage.candidates_token_count
            )
            
            # JSON 파싱
            try:
                repaired = repair_json(content)
                return json.loads(repaired)
            except:
                if "```json" in content:
                    content = content.split("```json")[1].split("```")[0]
                elif "```" in content:
                    content = content.split("```")[1].split("```")[0]
                return json.loads(content.strip())
                
        except Exception as e:
            logger.error(f"Gemini analysis error: {e}")
            return {
                "target": None,
                "indication": None,
                "summary": "Analysis failed",
                "relevance_score": 0.0,
                "ai_reasoning": f"Error: {str(e)}"
            }
    
    async def is_duplicate(self, title: str) -> bool:
        """knowledge_base에 이미 존재하는 논문인지 확인 (제목 기준)"""
        try:
            existing = supabase.table("knowledge_base")\
                .select("id")\
                .eq("title", title[:500])\
                .execute()
            return bool(existing.data)
        except:
            return False
    
    async def save_to_knowledge_base(
        self, 
        article: Dict, 
        analysis: Dict,
        drug_id: Optional[str] = None
    ) -> bool:
        """
        knowledge_base 스키마에 맞게 저장
        - source_type: 'PubMed'
        - title: 논문 제목
        - content: 초록 전문
        - summary: AI 요약
        - relevance_score: ADC 관련성
        - ai_reasoning: 중요성 판단 근거
        - rag_status: 'completed'
        """
        try:
            # AI 분석 결과로 요약 구성
            summary_parts = []
            if analysis.get("target"):
                summary_parts.append(f"Target: {analysis['target']}")
            if analysis.get("indication"):
                summary_parts.append(f"Indication: {analysis['indication']}")
            summary_parts.append(f"PMID: {article.get('pmid', 'N/A')}")
            
            summary_text = " | ".join(summary_parts)
            if analysis.get("summary"):
                summary_text += f"\n{analysis['summary']}"
            
            new_kb = {
                "source_type": "PubMed",
                "title": article["title"][:500],
                "content": article.get("abstract", "No abstract available."),
                "summary": summary_text[:1000],
                "relevance_score": analysis.get("relevance_score", 0.0),
                "source_tier": 1,  # PubMed = Tier 1 (학술 논문)
                "ai_reasoning": analysis.get("ai_reasoning", ""),
                "rag_status": "completed"  # AI 분석 완료
            }
            
            supabase.table("knowledge_base").insert(new_kb).execute()
            return True
            
        except Exception as e:
            logger.error(f"Failed to save to knowledge_base: {e}")
            return False
    
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
                # 중복 체크
                if await self.is_duplicate(article["title"]):
                    self.skipped_count += 1
                    continue
                
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
        mode: 'incremental' (신규만) 또는 'full' (전체)
        """
        from app.api.scheduler import update_job_status, is_cancelled
        
        if job_id:
            await update_job_status(job_id, status="running")
        
        logger.info(f"🚀 [PubMed Knowledge] Starting batch (size: {batch_size}, mode: {mode})")
        
        try:
            # 1. 대상 약물 가져오기
            limit = batch_size if mode == "incremental" else None
            drugs = await self.get_target_drugs(limit)
            
            if not drugs:
                logger.info("No target drugs found.")
                if job_id:
                    await update_job_status(job_id, status="completed", message="No drugs to process")
                return {"status": "completed", "processed": 0}
            
            if job_id:
                await update_job_status(job_id, records_found=len(drugs))
            
            # 2. 각 약물 처리
            total_saved = 0
            for idx, drug in enumerate(drugs):
                # 중단 체크
                if job_id and await is_cancelled(job_id):
                    logger.info("Job cancelled by user")
                    await update_job_status(job_id, status="stopped")
                    break
                
                saved = await self.process_single_drug(drug, job_id)
                total_saved += saved
                
                # 진행률 업데이트 (10개마다)
                if job_id and (idx + 1) % 10 == 0:
                    await update_job_status(
                        job_id, 
                        records_drafted=self.processed_count,
                        message=f"Processing {idx + 1}/{len(drugs)} drugs..."
                    )
                    logger.info(f"📊 Progress: {idx + 1}/{len(drugs)} drugs, {self.processed_count} papers saved")
                
                # Rate limiting between drugs
                await asyncio.sleep(self.REQUEST_DELAY)
            
            # 3. 완료
            result = {
                "status": "completed",
                "total_drugs": len(drugs),
                "papers_saved": self.processed_count,
                "skipped_duplicates": self.skipped_count,
                "errors": self.error_count,
                "failed_drugs": self.failed_drugs[:50]  # 최대 50개만 저장
            }
            
            if job_id:
                await update_job_status(
                    job_id,
                    status="completed",
                    records_drafted=self.processed_count,
                    completed_at=datetime.utcnow().isoformat(),
                    errors=self.failed_drugs[:20] if self.failed_drugs else [],
                    message=f"Saved {self.processed_count} papers from {len(drugs)} drugs"
                )
            
            logger.info(f"🎉 [PubMed Knowledge] Complete! Saved: {self.processed_count}, Errors: {self.error_count}")
            return result
            
        except Exception as e:
            logger.error(f"Batch processing error: {e}")
            if job_id:
                await update_job_status(job_id, status="failed", errors=[str(e)])
            return {"status": "failed", "error": str(e)}


# 싱글톤 인스턴스
pubmed_knowledge_service = PubMedKnowledgeService()
