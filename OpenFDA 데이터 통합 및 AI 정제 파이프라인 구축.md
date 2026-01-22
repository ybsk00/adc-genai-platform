[작업지시서] OpenFDA 데이터 통합 및 AI 정제 파이프라인 구축
1. 개요
목표: ClinicalTrials.gov와 동일한 아키텍처를 활용하여 OpenFDA 약물 승인 및 라벨 데이터를 golden_set_library 테이블로 통합.

모델: Gemini 2.5 Flash (SMILES Fallback 및 Chemist Persona 포함).

핵심 원칙: 기존 ClinicalTrials.gov 워커와 최대한 유사한 로직을 적용하여 개발 리소스 최소화.

2. 데이터 수집 (Ingestion) 로직
Endpoint: OpenFDA Drug Label API (또는 NDC/Adverse Events API).

수집 모드:

Full Import: ADC 관련 키워드(Antibody-drug conjugate 등)로 검색된 전체 데이터를 백그라운드에서 벌크 주입.

Daily Import: 신규 업데이트된 승인/라벨 정보만 한국 시간(KST) 새벽 4시에 정기 수집.

중복 처리: 약물 이름과 제조사를 기준으로 기존 레코드가 있다면 UPDATE(Upsert), 없으면 INSERT.

3. AI Refiner 고도화 (OpenFDA 맞춤형)
프롬프트 튜닝: 기존 임상용 프롬프트를 OpenFDA용으로 살짝 수정.

추출 항목: 승인 상태(Approval Status), 블랙박스 경고(Boxed Warning), 제조사 정보, 주요 부작용.

화학 구조 매핑:

Fuzzy Matching: 약물 이름으로 PubChem API 조회.

Gemini Fallback: 조회 실패 시 Gemini 2.5 Flash가 직접 SMILES 코드 생성 및 RDKit 검증.

출처 기록: enrichment_source 컬럼에 open_fda_api 또는 AI-Generated 명시.

4. UI/UX 요구사항 (Admin Dashboard)
AIRefinerStatusCard.tsx 및 관련 탭에 아래 기능 추가:

모드별 Run Now:

Run Batch (50): OpenFDA 미처리 데이터 중 50개만 즉시 분석.

Run Full (All): OpenFDA 모든 Pending 데이터를 전량 분석.

Run Daily Import: OpenFDA 최신 데이터를 즉시 수집(Trigger).

Spot Check: OpenFDA 데이터만 별도로 샘플링하여 AI 분석 결과를 눈으로 확인할 수 있는 기능 유지.

5. DB 스키마 유지 및 확장
동일 테이블 사용: golden_set_library.

데이터 격리: properties (JSONB) 필드 내에 fda_label_data 객체를 추가하여 기존 임상 데이터와 공존하도록 처리.

💡 사장님을 위한 기대 효과
이렇게 작업이 완료되면, 사장님은 하나의 화면에서 **"이 약물은 지금 임상이 어디까지 왔고(ClinicalTrials), FDA에서는 어떤 부작용 위험을 경고했는지(OpenFDA)"**를 한눈에 볼 수 있게 됩니다.

**"ClinicalTrials랑 똑같이, 대신 OpenFDA 데이터로 한 번 더 돌려놔!"**라고 이 지시서를 전달해 주십시오.

이제 사장님의 플랫폼은 단순한 정보 수집기를 넘어 **'글로벌 ADC 통합 인텔리전스'**로 진화하기 직전입니다. 개발자가 OpenFDA 워커를 세팅하는 동안, 우리는 이 데이터들을 활용해 어떤 리포트 기능을 먼저 만들지 논의해 볼까요?


한 번에 다 긁어오기 위한 "Full-Scan" 쿼리 (Backend용)
개발자가 API 엔드포인트에서 호출할 때 사용할 SQL 로직입니다.

SQL
-- 1. 아직 정제되지 않은(Pending) 모든 OpenFDA 데이터 조회
SELECT * FROM golden_set_library 
WHERE enrichment_source = 'open_fda_api' 
  AND rag_status = 'pending'
-- LIMIT 50;  <-- 이 부분을 제거하거나 5000으로 수정하도록 지시하십시오.
ORDER BY created_at ASC;
2. 🛠️ 대량 데이터 삽입 시 성능 최적화 (Upsert 쿼리)
OpenFDA API에서 수만 건을 가져올 때 DB가 뻗지 않게 하는 배치 삽입(Bulk Insert) 방식입니다.

SQL
-- 데이터 주입 시 중복은 업데이트하고 신규는 추가하는 로직
INSERT INTO golden_set_library (name, category, properties, enrichment_source, rag_status)
VALUES ($1, $2, $3, 'open_fda_api', 'pending')
ON CONFLICT (name, properties->>'manufacturer') -- 약물명과 제조사 기준 중복 체크
DO UPDATE SET 
    properties = golden_set_library.properties || EXCLUDED.properties, -- 기존 정보에 FDA 정보 병합
    updated_at = NOW();
3. 🚦 "Full" 모드 가동을 위한 개발자 지침
사장님이 원하시는 **"한번에 끝내기"**를 위해 개발자에게 아래 3가지를 확답받으세요.

페이지네이션 해제: OpenFDA API는 한 번 호출에 최대 99건 정도만 줍니다. 따라서 **"전체 데이터를 다 가져올 때까지 Loop(반복문)를 돌려라"**고 지시해야 합니다.

세마포어 확장: Gemini 2.5 Flash는 속도가 빠르므로, 동시 처리(Concurrency) 제한을 50개 이상으로 높여서 병렬로 정제하게 하세요.

서버 타임아웃 연장: 2,000건 이상을 한 번에 처리하려면 연결 시간이 길어질 수 있으므로, API 타임아웃 설정을 10분 이상으로 넉넉히 잡으라고 하십시오.

📊 조치 후 기대 지표
수집 속도: 50개씩 끊어서 클릭할 필요 없이, 단 한 번의 "Run Full" 클릭으로 전체 리스트가 채워집니다.

데이터 일관성: ClinicalTrials.gov와 OpenFDA 데이터가 SMILES를 중심으로 예쁘게 정렬됩니다.

**"50개 리밋은 테스트용일 뿐이다. 이제 실전이니 SQL 리밋 풀고, 오픈FDA 데이터 전량을 한 번에 밀어 넣어라!"**라고 말씀하시면 됩니다.

이제 OpenFDA까지 들어오면 정말 ADC 전문 포털로서의 면모를 갖추게 됩니다. 이 방대한 작업이 끝나는 대로, **투자자들에게 보여줄 [시장 선점 전략 리포트]**를 바로 뽑아볼까요?


중복 없이 전량 수집하기 위한 개발자 필독 로직
1. "건너뛰기(Skip)" 로직 도입
OpenFDA API는 skip 파라미터를 지원합니다. 이를 활용해 페이지를 넘기듯 데이터를 가져와야 합니다.

지시 사항: "한 번에 100개씩 가져오되, 다음 호출 때는 skip=100, 그다음은 skip=200 식으로 파라미터를 증가시키며 **전체 데이터(Total count)를 다 읽을 때까지 반복(Loop)**해라."

2. DB 레벨의 "중복 방지 방어막(Upsert)" 설치
API에서 중복된 데이터를 보내더라도 DB에 들어갈 때 걸러내야 합니다.

지시 사항: "golden_set_library 테이블에 데이터를 넣을 때 ON CONFLICT (name, enrichment_source) DO NOTHING 옵션을 써라. 이미 이름이 같은 OpenFDA 데이터가 있으면 아예 무시하고 넘어가서 DB 오염을 막아라."

3. "Daily"와 "Full"의 쿼리 차별화
Full 모드: skip을 계속 늘려가며 수천 건의 과거 승인 데이터를 전량 수집.

Daily 모드: OpenFDA의 last_updated 필드를 필터로 사용하여, **"어제 날짜 이후로 업데이트된 것"**만 콕 집어서 가져오도록 로직을 분리해라.

📝 개발자에게 전달할 [수정된 SQL 및 API 로직]
[OpenFDA 전량 수집 핵심 가이드]

Loop 수집: total 개수를 먼저 파악하고, limit=100과 skip 값을 루프 돌려 전수 조사할 것.

Filter 적용: search=openfda.substance_name:"antibody-drug conjugate" 처럼 키워드 필터를 API 레벨에서 걸어서 불필요한 일반 약물 데이터 소모를 방지할 것.

Upsert 필수: 동일 약물이 여러 번 수집되어도 golden_set_library에는 유니크하게 한 건만 남도록 처리할 것.