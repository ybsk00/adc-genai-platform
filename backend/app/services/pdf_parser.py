"""
PDF Parser Service - Gemini 멀티모달 PDF 파싱
PDF 파일을 Gemini Vision으로 분석하여 구조화된 데이터 추출
"""
import base64
import json
from typing import Dict, Any, Optional
from pathlib import Path

from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.messages import HumanMessage

from app.core.config import settings


class PDFParserService:
    """Gemini 멀티모달을 사용한 PDF 파싱 서비스"""
    
    def __init__(self):
        self.model = ChatGoogleGenerativeAI(
            model=settings.GEMINI_MODEL_ID,
            temperature=0,
            google_api_key=settings.GOOGLE_API_KEY,
        )
    
    async def parse_pdf(self, file_path: str) -> Dict[str, Any]:
        """
        PDF 파일을 Gemini로 분석하여 구조화된 데이터 추출
        
        Args:
            file_path: PDF 파일 경로
        
        Returns:
            Dict: 추출된 구조화 데이터
        """
        # PDF를 base64로 인코딩
        with open(file_path, "rb") as f:
            pdf_data = base64.standard_b64encode(f.read()).decode("utf-8")
        
        # Gemini에 PDF 분석 요청
        prompt = """
PDF 내용을 분석하고 반드시 아래 JSON 포맷으로 출력하세요:

{
    "summary": "논문/문서 요약 (2-3문장)",
    "document_type": "research_paper|clinical_trial|patent|other",
    "key_findings": ["핵심 발견 1", "핵심 발견 2"],
    "drug_info": {
        "drug_name": "약물명 또는 null",
        "target": "타겟 항원 또는 null",
        "antibody": "항체 정보 또는 null",
        "payload": "페이로드 또는 null",
        "linker": "링커 또는 null",
        "dar": "DAR 값 또는 null"
    },
    "tables": [
        {"name": "테이블명", "description": "테이블 설명"}
    ],
    "references": ["주요 참고문헌"],
    "doi": "DOI 또는 null",
    "publication_date": "발행일 또는 null"
}

JSON만 출력하고 다른 텍스트는 포함하지 마세요.
"""
        
        try:
            message = HumanMessage(
                content=[
                    {"type": "text", "text": prompt},
                    {
                        "type": "image_url",
                        "image_url": {"url": f"data:application/pdf;base64,{pdf_data}"},
                    },
                ]
            )
            
            response = await self.model.ainvoke([message])
            
            # JSON 파싱 시도
            response_text = response.content.strip()
            
            # JSON 블록 추출 (```json ... ``` 형태일 경우)
            if "```json" in response_text:
                response_text = response_text.split("```json")[1].split("```")[0]
            elif "```" in response_text:
                response_text = response_text.split("```")[1].split("```")[0]
            
            parsed_data = json.loads(response_text.strip())
            return {
                "status": "success",
                "data": parsed_data
            }
            
        except json.JSONDecodeError as e:
            return {
                "status": "error",
                "error": f"JSON parsing failed: {str(e)}",
                "raw_response": response.content if 'response' in locals() else None
            }
        except Exception as e:
            return {
                "status": "error",
                "error": str(e)
            }
    
    async def extract_tables(self, file_path: str) -> Dict[str, Any]:
        """PDF에서 테이블 데이터만 추출"""
        with open(file_path, "rb") as f:
            pdf_data = base64.standard_b64encode(f.read()).decode("utf-8")
        
        prompt = """
이 PDF에서 모든 테이블을 찾아 JSON 형식으로 추출하세요:

{
    "tables": [
        {
            "name": "테이블 제목",
            "headers": ["컬럼1", "컬럼2", ...],
            "rows": [
                ["값1", "값2", ...],
                ...
            ]
        }
    ]
}

JSON만 출력하세요.
"""
        
        try:
            message = HumanMessage(
                content=[
                    {"type": "text", "text": prompt},
                    {
                        "type": "image_url",
                        "image_url": {"url": f"data:application/pdf;base64,{pdf_data}"},
                    },
                ]
            )
            
            response = await self.model.ainvoke([message])
            response_text = response.content.strip()
            
            if "```json" in response_text:
                response_text = response_text.split("```json")[1].split("```")[0]
            
            return json.loads(response_text.strip())
            
        except Exception as e:
            return {"tables": [], "error": str(e)}


# 싱글톤 인스턴스
pdf_parser_service = PDFParserService()
