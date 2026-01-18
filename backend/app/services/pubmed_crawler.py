"""
PubMed Crawler Service - PubMed 자동 수집 봇
ADC 관련 논문을 자동으로 수집하고 AI Gatekeeper로 필터링
"""
import httpx
from typing import List, Dict, Any, Optional
from datetime import datetime, timedelta
import asyncio

from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.prompts import ChatPromptTemplate

from app.core.config import settings


class PubMedCrawlerService:
    """PubMed 논문 자동 수집 서비스"""
    
    PUBMED_API_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    DEFAULT_QUERY = '("antibody drug conjugate" OR "ADC") AND ("linker" OR "payload")'
    
    def __init__(self):
        self.gemini = ChatGoogleGenerativeAI(
            model=settings.GEMINI_MODEL_ID,
            temperature=0,
            google_api_key=settings.GOOGLE_API_KEY,
        )
    
    async def search_pubmed(
        self, 
        query: str = None, 
        max_results: int = 20,
        days_back: int = 7
    ) -> List[str]:
        """
        PubMed에서 최근 논문 검색
        
        Args:
            query: 검색 쿼리 (기본값: ADC 관련)
            max_results: 최대 결과 수
            days_back: 며칠 전까지 검색할지
        
        Returns:
            List[str]: PubMed ID 목록
        """
        query = query or self.DEFAULT_QUERY
        
        # 날짜 범위 설정
        end_date = datetime.now()
        start_date = end_date - timedelta(days=days_back)
        
        search_url = f"{self.PUBMED_API_BASE}/esearch.fcgi"
        params = {
            "db": "pubmed",
            "term": query,
            "retmax": max_results,
            "retmode": "json",
            "datetype": "pdat",
            "mindate": start_date.strftime("%Y/%m/%d"),
            "maxdate": end_date.strftime("%Y/%m/%d"),
            "sort": "date"
        }
        
        async with httpx.AsyncClient() as client:
            response = await client.get(search_url, params=params)
            data = response.json()
        
        return data.get("esearchresult", {}).get("idlist", [])
    
    async def fetch_article_details(self, pmid: str) -> Dict[str, Any]:
        """
        PubMed 논문 상세 정보 조회
        
        Args:
            pmid: PubMed ID
        
        Returns:
            Dict: 논문 상세 정보
        """
        fetch_url = f"{self.PUBMED_API_BASE}/efetch.fcgi"
        params = {
            "db": "pubmed",
            "id": pmid,
            "retmode": "xml",
            "rettype": "abstract"
        }
        
        async with httpx.AsyncClient() as client:
            response = await client.get(fetch_url, params=params)
            xml_content = response.text
        
        # XML 파싱 (간단한 추출)
        import re
        
        title_match = re.search(r"<ArticleTitle>(.+?)</ArticleTitle>", xml_content, re.DOTALL)
        abstract_match = re.search(r"<AbstractText.*?>(.+?)</AbstractText>", xml_content, re.DOTALL)
        doi_match = re.search(r'<ArticleId IdType="doi">(.+?)</ArticleId>', xml_content)
        
        return {
            "pmid": pmid,
            "title": title_match.group(1).strip() if title_match else None,
            "abstract": abstract_match.group(1).strip() if abstract_match else None,
            "doi": doi_match.group(1) if doi_match else None,
            "source": "pubmed",
            "source_id": pmid
        }
    
    async def ai_gatekeeper(self, article: Dict[str, Any]) -> Dict[str, Any]:
        """
        AI Gatekeeper - 논문이 ADC 관련인지 검수
        
        Args:
            article: 논문 정보 (title, abstract 포함)
        
        Returns:
            Dict: 검수 결과
        """
        if not article.get("abstract"):
            return {"is_relevant": False, "reason": "No abstract available"}
        
        prompt = ChatPromptTemplate.from_messages([
            ("system", """You are an AI gatekeeper that evaluates if a scientific article is relevant to ADC (Antibody-Drug Conjugate) research.

Evaluate based on these criteria:
1. Is the article about ADC drugs or ADC development?
2. Does it discuss linkers, payloads, or conjugation chemistry?
3. Is it relevant to ADC toxicity, efficacy, or clinical trials?

Respond in JSON format:
{
    "is_relevant": true/false,
    "confidence": 0.0-1.0,
    "reason": "Brief explanation",
    "category": "structure|toxicity|clinical|efficacy|other|not_relevant"
}
"""),
            ("user", """Title: {title}

Abstract: {abstract}

Is this article relevant to ADC research?""")
        ])
        
        chain = prompt | self.gemini
        
        try:
            response = await chain.ainvoke({
                "title": article.get("title", ""),
                "abstract": article.get("abstract", "")
            })
            
            import json
            response_text = response.content.strip()
            if "```json" in response_text:
                response_text = response_text.split("```json")[1].split("```")[0]
            
            result = json.loads(response_text.strip())
            result["pmid"] = article.get("pmid")
            return result
            
        except Exception as e:
            return {
                "is_relevant": False,
                "reason": f"Evaluation error: {str(e)}",
                "pmid": article.get("pmid")
            }
    
    async def crawl_and_filter(
        self, 
        query: str = None, 
        max_results: int = 20,
        days_back: int = 7
    ) -> Dict[str, Any]:
        """
        논문 수집 + AI 필터링 전체 프로세스
        
        Returns:
            Dict: 수집 및 필터링 결과
        """
        # 1. PubMed 검색
        pmids = await self.search_pubmed(query, max_results, days_back)
        
        if not pmids:
            return {
                "status": "success",
                "total_found": 0,
                "relevant": [],
                "filtered_out": []
            }
        
        # 2. 상세 정보 조회
        articles = []
        for pmid in pmids:
            try:
                article = await self.fetch_article_details(pmid)
                articles.append(article)
                await asyncio.sleep(0.5)  # Rate limiting
            except Exception as e:
                print(f"Error fetching {pmid}: {e}")
        
        # 3. AI Gatekeeper 필터링
        relevant = []
        filtered_out = []
        
        for article in articles:
            result = await self.ai_gatekeeper(article)
            
            if result.get("is_relevant"):
                relevant.append({
                    **article,
                    "gatekeeper_result": result
                })
            else:
                filtered_out.append({
                    "pmid": article.get("pmid"),
                    "title": article.get("title"),
                    "reason": result.get("reason")
                })
        
        return {
            "status": "success",
            "total_found": len(pmids),
            "total_processed": len(articles),
            "relevant_count": len(relevant),
            "filtered_out_count": len(filtered_out),
            "relevant": relevant,
            "filtered_out": filtered_out
        }


# 싱글톤 인스턴스
pubmed_crawler = PubMedCrawlerService()
