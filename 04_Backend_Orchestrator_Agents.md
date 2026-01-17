04. Backend Design: Orchestrator & AI Agents
Document ID: BE-01 Role: The "Brain" of the Platform (Reasoning & Processing) Tech Stack: Python 3.10+, FastAPI, LangChain / LangGraph, Pydantic External APIs: Google Vertex AI (AlphaFold), NVIDIA BioNeMo, Perplexity (Sonar), OpenAI (GPT-4o)

1. Orchestration Workflow (업무 흐름도)
우리는 속도를 위해 **병렬 처리(Parallel Execution)**와 **순차 처리(Sequential)**를 혼합합니다.

코드 스니펫

graph TD
    Start([User Input: Sequence + SMILES]) --> Router{Orchestrator}
    
    subgraph "Phase 1: Parallel Analysis (동시 작업)"
        Router --> A1[Agent 1: Structure]
        Router --> A2[Agent 2: Toxicology]
        Router --> A3[Agent 3: Patent]
        Router --> A4[Agent 4: Competitor]
    end
    
    A2 --> A5
    A4 --> A5
    
    subgraph "Phase 2: Synthesis (종합)"
        A5[Agent 5: Clinical Design]
        A1 & A3 & A5 --> A6[Agent 6: Report Writer]
    end
    
    A6 --> PDF_Gen[PDF Generator] --> End([Final Report])
Logic: '독성(A2)'과 '경쟁사(A4)' 분석이 끝나야만, 그걸 바탕으로 '임상 프로토콜(A5)'을 짤 수 있습니다. (Dependency)

State Management: 모든 에이전트는 ADCState라는 하나의 공유 메모리 객체를 읽고 씁니다.

2. State Schema (데이터 구조)
에이전트끼리 주고받는 **작업 지시서(Context)**의 형태입니다. (Pydantic Model)

Python

class ADCState(BaseModel):
    # Input
    job_id: str
    antibody_seq: str
    payload_smiles: str
    target_name: str  # e.g., "LIV-1"
    
    # Analysis Results (채워질 내용들)
    structure_data: dict = Field(default={})  # 3D coords, PDB path
    toxicity_data: dict = Field(default={})   # LogP, risk_flags
    patent_data: dict = Field(default={})     # FTO status
    competitor_data: dict = Field(default={}) # Competitor list
    clinical_data: dict = Field(default={})   # Draft protocol
    
    # Final Output
    final_report_json: dict = Field(default={})
    status: str = "processing" # processing, completed, failed
3. The 6-Agent Squad Specifications
각 에이전트의 페르소나(Persona), 사용 도구(Tools), **시스템 프롬프트(System Prompt)**입니다.

🕵️‍♂️ Agent 1. Bio-Structure Agent (구조 분석가)
Role: 항체와 약물이 물리적으로 잘 붙는지, 엉키지는 않는지 3D 구조 관점에서 분석.

Tools:

Google Vertex AI (AlphaFold): 항체 3D 구조 예측.

RDKit: 페이로드(약물) 3D 구조 생성.

System Prompt:

Plaintext

You are an expert Structural Biologist specializing in ADCs.
Your task is to analyze the 3D compatibility between the provided antibody sequence and the payload.
1. Check for steric hindrance in the CDR regions.
2. Identify surface-exposed cysteine residues for conjugation.
3. Output the path to the generated .pdb file and a summary of structural stability.
☣️ Agent 2. Toxicology Agent (독성학자)
Role: 약물의 화학적 특성을 보고 부작용 위험을 예측. (가장 중요)

Tools:

Vector Search (RAG): Golden Set DB에서 유사 약물의 실패 사례 검색.

RDKit: LogP(소수성), TPSA(표면적) 계산.

System Prompt:

Plaintext

You are a Senior Toxicologist.
Analyze the payload's SMILES code and predict potential risks.
- Calculate LogP. If LogP > 3.5, flag "Aggregation Risk".
- Search the Golden Set DB for similar payloads (e.g., MMAE, DXd).
- If similar payloads caused "Ocular Toxicity" or "Neutropenia" in trials, flag it as HIGH RISK.
⚖️ Agent 3. Patent Agent (변리사)
Role: 4대 플랫폼(엔허투, 시젠 등)의 특허 침해 여부 판단.

Tools: Tavily Search API (Google Search), Internal Patent DB.

System Prompt:

Plaintext

You are an IP Attorney specializing in Biotech.
Compare the user's linker-payload structure against major ADC patents (Daiichi Sankyo, Seagen).
- Identify if the linker sequence (e.g., GGFG) matches patented sequences.
- Provide a "Traffic Light" status: Green (Safe/Expired), Yellow (Modify), Red (Infringement).
🔭 Agent 4. Competitor Agent (경쟁사 분석가)
Role: 해당 타겟(LIV-1)을 연구하는 다른 회사를 찾고 비교.

Tools: Perplexity API (sonar-medium-online) - 실시간 웹 검색 최강.

System Prompt:

Plaintext

You are a Market Intelligence Analyst.
Search for "ADC drugs targeting {target_name}" in clinical trials and news.
- List top 3 competitors.
- Compare their development stage (Phase 1/2/3) and payloads.
- Identify their weaknesses (e.g., "High toxicity reported in Phase 1").
👩‍⚕️ Agent 5. Clinical Design Agent (임상 기획자)
Role: 독성 데이터(A2)와 경쟁사 데이터(A4)를 보고 임상 1상 설계를 제안.

Inputs: toxicity_data, competitor_data from previous agents.

System Prompt:

Plaintext

You are a Clinical Development Director.
Based on the toxicology report (Agent 2) and competitor landscape (Agent 4), design a draft Phase 1 protocol.
- Suggest 'Starting Dose' (e.g., if tox is high, start low at 0.5 mg/kg).
- Define 'Inclusion Criteria' (e.g., HER2-low patients).
✍️ Agent 6. Report Writer (편집장)
Role: 모든 에이전트의 JSON 결과를 취합하여, 최종 PDF용 JSON을 생성하고 등급(Grade)을 매김.

LLM: GPT-4o (Reasoning 능력이 가장 좋음).

System Prompt:

Plaintext

You are the Chief Scientific Officer (CSO).
Synthesize the reports from all 5 agents into a coherent "Investment Memo".
- Assign a final grade (S/A/B/C/F).
- Write a "Verdict" (Go / No-Go).
- Highlight the biggest risk and the best opportunity.
- Output strictly in JSON format matching the Report Template schema.
4. Error Handling & Fallback (에러 처리)
AI 에이전트 하나가 죽어도 전체 리포트가 실패하면 안 됩니다.

Retry Policy: 각 에이전트는 최대 3회 재시도 (Exponential Backoff).

Graceful Degradation:

만약 Agent 4 (Competitor)가 Perplexity 오류로 실패하면?

전체를 멈추지 않고, 경쟁사 섹션에 **"Data currently unavailable (External API Error)"**라고 적고 리포트를 발행합니다.

단, Agent 2 (Toxicology)가 실패하면 Critical Error로 간주하고 작업을 중단합니다. (독성 없는 리포트는 쓰레기니까요.)

5. Development Checklist
[ ] LangGraph 설치 및 노드(Node)/엣지(Edge) 구성.

[ ] Google Cloud Run에 배포할 Dockerfile에 BioNeMo, RDKit 의존성 추가.

[ ] Environment Variables 설정:

OPENAI_API_KEY

PERPLEXITY_API_KEY

BIONEMO_API_KEY

TAVILY_API_KEY

[ ] Mock Test: 실제 API를 부르지 않고 가짜 데이터(Mock)로 에이전트 흐름이 잘 넘어가는지 단위 테스트 작성.