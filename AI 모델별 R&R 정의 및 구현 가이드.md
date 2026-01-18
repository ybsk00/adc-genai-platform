[작업지시서] AI 모델별 R&R 정의 및 구현 가이드
문서 번호: SPEC-AI-01 작성 일자: 2026-01-18 수신: 개발팀 목표: 비용 효율화 및 처리 속도 향상을 위한 Multi-LLM Hybrid Architecture 구현

1. 개요 (Overview)
기존 단일 모델(GPT) 의존도에서 벗어나, 작업의 성격에 따라 가장 적합한 AI 모델을 골라 쓰는 라우팅(Routing) 시스템을 구축합니다.

Gemini 2.0 Flash: 막대한 양의 논문/특허 데이터 처리, 단순 요약, 속도가 중요한 RAG.

GPT-4o: 고도의 논리적 추론, 최종 의사결정(Verdict), 정교한 JSON 출력.

Perplexity: 실시간 최신 뉴스, 시장 동향 파악 (Web Search).

2. 모델별 상세 역할 (Detailed Roles)
🤖 1. Google Gemini 2.0 Flash
별명: The Fast Reader (속독왕)

역할 정의:

RAG (검색 증강 생성): 벡터 DB에서 찾아낸 수십 페이지의 논문 조각(Chunk)을 읽고 핵심 정보를 추출.

초기 데이터 파싱: PDF 파일에서 표(Table) 데이터나 수치를 뽑아내는 단순 노동.

1~4번 에이전트(실무자): 구조 분석, 독성 수치 확인 등 '팩트'를 찾는 에이전트들의 두뇌.

선정 이유:

압도적인 Context Window (긴 글을 한 번에 읽음).

GPT-4o 대비 속도가 매우 빠르고 비용이 저렴함.

적용 코드 위치: agent_structure.py, agent_toxicology.py

🧠 2. OpenAI GPT-4o
별명: The Logic Master (똑똑한 팀장)

역할 정의:

Agent 6 (리포트 작성자): 각 에이전트가 조사한 내용을 종합하여 최종 보고서 작성.

오케스트레이터: "독성 정보가 부족하니 다시 찾아와" 같은 판단 및 지시.

JSON Enforcement: 프론트엔드에서 깨지지 않는 완벽한 JSON 규격 생성.

선정 이유:

복잡한 지시사항 이행 능력(Instruction Following)이 세계 최고 수준.

최종 결과물의 문장 퀄리티가 가장 자연스러움.

적용 코드 위치: agent_orchestrator.py, agent_report_writer.py

🔎 3. Perplexity (Sonar)
별명: The News Hunter (기자)

역할 정의:

외부 검색: "최근 화이자의 ADC 인수 뉴스 찾아줘"와 같은 실시간 정보 수집.

시장 동향 파악: 학습 데이터에 없는 최신 트렌드 조사.

적용 코드 위치: tool_market_search.py (Worker)

3. 구현 가이드 (Implementation Spec)
개발자는 LangChain을 사용하여 아래와 같이 모델을 스위칭해야 합니다.

3.1 환경변수 설정 (backend/.env)
Properties

# Google (Gemini 2.0 Flash) - 분석용
GOOGLE_API_KEY="AIzaSy..."
GEMINI_MODEL_ID="gemini-2.0-flash-exp"

# OpenAI (GPT-4o) - 종합용
OPENAI_API_KEY="sk-proj-..."
OPENAI_MODEL_ID="gpt-4o"

# Perplexity - 검색용
PERPLEXITY_API_KEY="pplx-..."
3.2 Python 코드 예시 (Pseudo-code)
Case A: 논문 분석 에이전트 (Gemini 사용)

Python

from langchain_google_genai import ChatGoogleGenerativeAI

# 분석가는 속도와 긴 문맥이 중요 -> Gemini 2.0 Flash 호출
def get_toxicology_analysis(context_text):
    llm = ChatGoogleGenerativeAI(
        model="gemini-2.0-flash-exp",
        temperature=0.1, # 팩트 위주라 창의성 낮춤
        google_api_key=os.getenv("GOOGLE_API_KEY")
    )
    
    prompt = f"다음 논문 내용을 바탕으로 IC50 값을 추출해: {context_text}"
    return llm.invoke(prompt)
Case B: 최종 리포트 작성 에이전트 (GPT-4o 사용)

Python

from langchain_openai import ChatOpenAI

# 최종 작성자는 논리력과 포맷 준수가 중요 -> GPT-4o 호출
def generate_final_report(all_agent_outputs):
    llm = ChatOpenAI(
        model="gpt-4o",
        temperature=0.5, # 자연스러운 문장을 위해 약간 높임
        openai_api_key=os.getenv("OPENAI_API_KEY")
    )
    
    prompt = f"""
    아래 전문가들의 분석을 종합해서 CEO 보고용 리포트를 작성해.
    반드시 JSON 포맷을 지켜야 해.
    
    [분석 데이터]
    {all_agent_outputs}
    """
    return llm.invoke(prompt)
4. 핵심 요청 사항 (Summary)
비용 절감: 텍스트 양이 많은(Heavy) 단순 분석 작업은 무조건 Gemini 2.0 Flash로 돌려주세요.

품질 보장: 고객이 보는 최종 결과물(Final Output)은 비싸더라도 GPT-4o를 써주세요.

검색: 실시간 정보가 필요할 때만 Perplexity API를 호출하세요.