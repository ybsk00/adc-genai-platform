사용자님의 ADC 플랫폼이 "똑똑한 컨설턴트"가 되려면, **Golden Set(논문, 특허, 임상 데이터)**을 AI가 이해하기 가장 좋은 형태로 잘게 썰어서(Chunking) 먹여줘야 합니다.

07_RAG_Pipeline_Chunking.md 문서는 바이오/화학 데이터의 특성을 고려한 특수 파싱 및 청킹 전략을 정의합니다. 일반적인 텍스트 청킹과는 차원이 다른, 전문가 수준의 데이터 처리 설계입니다.

07. RAG Pipeline: Parsing & Semantic Chunking
Document ID: DATA-01 Role: ETL (Extract, Transform, Load) for Unstructured Data Tech Stack: Python, LlamaParse (PDF Parser), LangChain, OpenAI Embeddings (text-embedding-3-small)

1. Pipeline Overview (데이터 흐름)
단순히 텍스트만 추출하는 것이 아니라, **표(Table)**와 **화학 구조(Chemicals)**를 별도로 처리하는 것이 핵심입니다.

코드 스니펫

graph TD
    Input[Source Documents: PDF/HTML] --> Router{Document Type?}
    
    subgraph "PDF Parsing (LlamaParse)"
        Router -->|Research Paper/Patent| P1[Extract Text]
        Router -->|Tables & Figures| P2[Extract Tables to Markdown]
    end
    
    subgraph "Semantic Chunking"
        P1 --> C1[Text Splitter (Recursive)]
        P2 --> C2[Table Summarizer (LLM)]
    end
    
    subgraph "Enrichment"
        C1 & C2 --> E1[Metadata Tagging]
        E1 --> E2[Chemical Entity Recognition (NER)]
    end
    
    E2 --> Embed[Vector Embedding] --> DB[(Supabase pgvector)]
2. Parsing Strategy (파싱 전략)
바이오 논문은 **"복잡한 표"**와 **"화학식"**이 생명입니다. 일반적인 파서(PyPDF2 등)는 이걸 다 깨먹습니다.

2.1 PDF Parser Selection: LlamaParse
Why: 표(Table)를 망가뜨리지 않고 Markdown 포맷으로 깔끔하게 변환해 주는 현재 최고의 파서입니다.

Settings:

parsing_instruction: "This is a biotech research paper. Preserve all tables related to IC50 values and clinical trial results in Markdown format."

2.2 Table Handling (표 처리 - 핵심)
논문 속 '독성 데이터 표'는 벡터 검색에서 검색이 잘 안 됩니다. 이를 해결하기 위해 LLM 요약을 거칩니다.

Raw Table: | Drug | IC50 (nM) | Status | | :--- | :--- | :--- | | MMAE | 0.5 | Ph2 |

Transformation (LLM 요약):

"The drug MMAE shows an IC50 value of 0.5 nM and is currently in Phase 2 clinical trials."

전략: 원본 표(Markdown)와 요약문(Text)을 함께 청킹하여 저장합니다.

3. Semantic Chunking Logic (청킹 로직)
단순히 500자씩 자르면 맥락이 끊깁니다. **"의미 단위"**로 잘라야 합니다.

3.1 Hierarchical Chunking (계층적 청킹)
Parent Chunk (큰 덩어리): 논문의 한 섹션 전체 (예: "3. Results").

용도: LLM에게 답변 생성 시 전체 맥락 제공.

Child Chunk (작은 덩어리): 구체적인 실험 결과 문단.

용도: 벡터 검색(RAG)의 정확도 향상.

3.2 Metadata Tagging (메타데이터 태깅)
청크만 저장하면 나중에 정밀 검색이 불가능합니다. 청크마다 태그를 붙입니다.

Drug Name: MMAE, Trastuzumab (NER로 자동 추출)

Document Type: Patent, Clinical Trial, Paper

Year: 2024

Target: LIV-1

💡 검색 예시: "LIV-1 타겟(Target)인 2023년 이후(Year) 특허(Type)만 찾아줘." (이게 가능해집니다.)

4. Embedding Model (임베딩 모델)
Model: OpenAI text-embedding-3-small

Dimension: 1536

Why: 가성비가 가장 좋고, 한국어와 영어가 섞여 있어도 성능이 준수합니다.

Optimization: Matryoshka Representation Learning 기법을 지원하므로, 나중에 차원을 줄여서 속도를 높일 수도 있습니다.

5. Development Specs (구현 상세)
개발자가 작성해야 할 파이썬 코드(data_pipeline/parser/chunker.py)의 로직입니다.

Python

# Pseudo-code for Semantic Chunking

def process_document(pdf_path):
    # 1. Parse PDF to Markdown
    markdown_text = llama_parse.parse(pdf_path)
    
    # 2. Split by Headers (# Introduction, # Methods...)
    sections = split_by_markdown_headers(markdown_text)
    
    chunks = []
    for section in sections:
        # 3. Extract Chemical Names (MMAE, DXd...) using Regex or NER
        chemicals = extract_chemical_entities(section.text)
        
        # 4. Create Chunk Object
        chunk = {
            "content": section.text,
            "metadata": {
                "source": pdf_path,
                "chemicals": chemicals, # ["MMAE", "LIV-1"]
                "section": section.header
            }
        }
        chunks.append(chunk)
        
    # 5. Embedding & Save to DB
    vector_store.add_documents(chunks)
6. Development Checklist
[ ] LlamaParse API Key 발급 및 연동 테스트.

[ ] Table Summarizer 프롬프트 작성 (표를 글로 풀어주는 LLM 함수).

[ ] Supabase pgvector 인덱스(hnsw) 생성 확인.

[ ] Chemical NER (약물 이름 추출) 간단한 리스트 매칭 로직 구현.

💡 작성자의 코멘트
이 문서는 **"AI가 논문을 읽을 때 표를 놓치지 않게 하라"**는 것이 핵심입니다. 특히 2.2 Table Handling은 경쟁사들이 잘 못하는 부분(대부분 텍스트만 긁음)이므로, 우리 플랫폼의 검색 정확도를 비약적으로 높여줄 무기가 될 것입니다.