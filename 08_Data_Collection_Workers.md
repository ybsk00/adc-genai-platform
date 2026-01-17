08. Data Collection Workers & AutomationDocument ID: DATA-02Role: Continuous Data Ingestion (The "Fuel" Supply)Tech Stack: Python, Celery (Scheduler), ClinicalTrials.gov API, Perplexity API, BioPython (PubMed)1. System Overview (수집 공장 구조도)데이터 수집은 **'한 번에 왕창 붓기(Initial)'**와 **'조금씩 자주 업데이트하기(Incremental)'**로 나뉩니다.코드 스니펫graph LR
    Scheduler[Celery Beat (Cron)] -->|Trigger| WorkerNode{Worker Swarm}
    
    subgraph "Sources"
        WorkerNode -->|API| S1[ClinicalTrials.gov]
        WorkerNode -->|Library| S2[PubMed / BioRxiv]
        WorkerNode -->|Search| S3[Perplexity (News/Patents)]
    end
    
    S1 & S2 & S3 -->|Raw JSON/PDF| Staging[Raw Data Storage (S3)]
    Staging -->|Event Trigger| Parser[RAG Pipeline (Doc 07)]
    Parser -->|Update| DB[(Golden Set DB)]
2. Strategy A: Initial Seeding (초기 데이터 붓기)서비스 런칭 전, 텅 빈 DB를 채우기 위한 1회성 대규모 작업입니다.2.1 Clinical Trials (임상 데이터)Target: ClinicalTrials.gov 
(미국 NIH)Method: Direct JSON Download (크롤링 필요 없음)이 사이트는 전체 데이터를 JSON으로 덤프 떠서 제공합니다.Filter:Keyword: 
"Antibody Drug Conjugate" OR "ADC"Status: Recruiting, Active, Completed, 
TerminatedAction: 약 2,000건의 ADC 관련 임상 데이터를 다운로드하여 DB의 golden_set_library 테이블에 적재.2.2 Approved Drugs (승인된 약물 족보)Target: FDA Approved ADCs List 
(엔허투, 캐드실라 등 15종+)Action: 이 15개 약물은 관리자가 수동으로 아주 정밀하게 입력하는 것을 권장합니다.이유: 우리 서비스의 기준점(Benchmark)이 되므로, 1%의 오류도 없어야 합니다.

3. Strategy B: Incremental Update (자동 업데이트 워커)서비스 런칭 후, 매일/매주 새로운 정보를 물어오는 봇(Bot)들입니다.🕵️‍♂️ Worker 1. The Clinical Watcher 
(임상 감시자)Frequency: Weekly (매주 월요일 새벽 3시)API: ClinicalTrials.gov API v2Logic:지난주 이후 LastUpdatePostDate가 변경된 
ADC 임상 검색.만약 "Status"가 Phase 1 → Phase 2로 바뀌었다?[Alert] 관리자에게 알림을 보내고, DB 업데이트.[Value] 사용자 리포트의 "경쟁사 현황"이 자동으로 최신화됨.🔬 
Worker 2. The Paper Hunter (논문 사냥꾼)Frequency: DailyAPI: BioPython.Entrez (PubMed API)Query: ("Antibody-Drug Conjugate"[Title/Abstract]) 
AND ("2026/01"[Date - Publication])Logic:어제 나온 ADC 관련 신규 논문 초록(Abstract) 수집.LLM이 "새로운 링커 기술인가?" 판단.맞다면 RAG 파이프라인(Doc 07)
으로 보내서 벡터 DB에 저장.📰 Worker 3. The Market Spy (뉴스/동향 스파이)Frequency: Daily (매일 아침 8시)API: Perplexity API (or SerpApi)Query: 
"Latest ADC biotech deals and acquisitions this week"Logic:"화이자(Pfizer)가 시젠(Seagen) 인수" 같은 빅뉴스 감지.이 내용은 DB 저장보다는 **대시보드의 '뉴스 피드'**에 즉시 노출.

4. Implementation Specs (코드 로직)개발자가 작성해야 할 worker/clinical_crawler.py의 핵심 로직입니다.Pythonimport requests
from app.db import supabase

def fetch_new_trials():
    url = "https://clinicaltrials.gov/api/v2/studies"
    params = {
        "query.term": "Antibody Drug Conjugate",
        "filter.lastUpdatePostDate": "gte:2026-01-10", # 지난주
        "format": "json"
    }
    
    response = requests.get(url, params=params)
    studies = response.json()
    
    for study in studies['studies']:
        # 데이터 정제 (Extract Fields)
        trial_id = study['protocolSection']['identificationModule']['nctId']
        phase = study['protocolSection']['designModule']['phases'][0]
        drug_name = extract_drug_name(study) # 별도 NLP 함수
        
        # DB Upsert (없으면 넣고, 있으면 업데이트)
        supabase.table('golden_set_library').upsert({
            "source": trial_id,
            "name": drug_name,
            "properties": {"phase": phase},
            "category": "clinical_trial"
        }).execute()
5. Scheduling (Crontab 설정)백엔드 서버(Cloud Run)가 아니라, 별도의 스케줄러 인스턴스나 Cloud Scheduler에서 돕니다.Worker NameSchedule (Cron)PriorityPurposemarket_spy0 8 * * * (매일 08:00)High대시보드 뉴스 피드 갱신paper_hunter0 2 * * * (매일 02:00)Medium신기술 RAG DB 적재clinical_watcher0 3 * * 1 (월요일 03:00)Low임상 단계 변경 추적6. Development Checklist[ ] ClinicalTrials.gov API 연동 테스트 (Rate Limit 확인).[ ] PubMed API 이메일 등록 (연구용 무료 사용을 위해 필수).[ ] Perplexity API Key 발급 및 크레딧 충전.[ ] Celery Beat 설정 (또는 Google Cloud Scheduler + Pub/Sub 구성 추천).💡 작성자의 코멘트이 08번 문서 덕분에 사용자님의 플랫폼은 **"가만히 있어도 똑똑해지는 시스템"**이 됩니다.초기 데이터(Initial): 처음에는 FDA 승인 약물 위주로 **소수정예(High Quality)**로 채우십시오.자동 수집(Incremental): 
욕심내서 모든 걸 긁으려 하지 말고, **"임상 단계 변경(Phase Change)"**만 잘 추적해도 연구원들에겐 엄청난 가치가 있습니다.