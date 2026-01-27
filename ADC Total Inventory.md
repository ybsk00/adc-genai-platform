[작업지시서] ADC Total Inventory 고도화 및 데이터 정제 자동화
1. 인프라 및 수집 전략 (Infra & VPN)
VPN 연결 및 로테이션:

ExpressVPN을 미국 동부(New York, New Jersey) 서버에 우선 연결해라.

수집량 500~1,000건 단위로 미국 내 타 도시(Chicago, Miami 등)로 IP 로테이션을 수행하여 차단을 회피해라.

차단 방어 로직:

HTTP 403(Forbidden) 에러 발생 시 즉시 수집을 중단하고 로그를 남긴 뒤, IP 교체 후 재시도하는 에러 핸들링을 강화해라.

Sleeping 타임은 현재의 안전 모드(12~20초)를 유지하되, VPN 안정성이 확인되면 점진적으로 최적화해라.

2. 데이터베이스 및 수집 필드 (Database)
메뉴 명칭 변경: 어드민 내 기존 'Seed' 메뉴를 **'Total Inventory'**로 변경해라.

테이블 통합 및 매핑:

Creative Biolabs(CBL)의 데이터를 commercial_reagents 테이블에 통합해라.

CRITICAL: ambeed_cat_no 필드에 CBL 제품번호(ADC-P-xxx)를 매핑하고, source_name을 'Creative Biolabs'로 명시해라.

cas_number, molecular_weight, classification 등 상세 제원을 누락 없이 파싱해라.

3. AI 기반 데이터 정제 자동화 (AI Refinement)
CAS 기반 SMILES 자동 채우기:

cas_number는 있으나 smiles_code가 없는 데이터를 추출하여 PubChem API를 통해 CanonicalSMILES를 자동으로 수집하는 스크립트를 가동해라.

AI 검증 프로세스:

API로 가져온 SMILES가 실제 제품명과 일치하는지 LLM(AI)이 2차 검증하게 하고, 검증 완료 시 ai_refined = true로 마킹해라.

오류 로그: 보정 실패 시 summary 필드에 사유를 기록하여 추후 수동 검토가 가능하게 해라.

4. 관리자 화면(Admin UI) 구성
인벤토리 대시보드 구축:

Tab 1 (Antibodies): antibody_library 리스트 (ID, Cat No, Target, Host Species 노출).

Tab 2 (Reagents): commercial_reagents 리스트 (Name, CAS, SMILES 유무, Source 노출).

시각적 가이드:

smiles_code가 비어 있는 행은 빨간색 또는 경고 아이콘으로 표시하여 사장님이 정제 현황을 직관적으로 파악하게 해라.

상세 보기 클릭 시 오른쪽 패널에서 전체 필드를 확인할 수 있는 Master-Detail View를 적용해라.


db스키마

create table public.antibody_library (
  id uuid not null default gen_random_uuid (),
  product_name text not null,
  cat_no text not null,
  antibody_format text null,
  host_species text null,
  isotype text null,
  related_disease text null,
  full_spec jsonb null,
  source_url text null,
  source_name text null default 'Creative Biolabs'::text,
  embedding public.vector null,
  summary text null,
  crawled_at timestamp with time zone null default now(),
  updated_at timestamp with time zone null default now(),
  constraint antibody_library_pkey primary key (id),
  constraint antibody_library_cat_no_key unique (cat_no)
) TABLESPACE pg_default;

create trigger handle_updated_at_antibody BEFORE
update on antibody_library for EACH row
execute FUNCTION moddatetime ('updated_at');



create table public.commercial_reagents (
  id uuid not null default gen_random_uuid (),
  ambeed_cat_no text not null,
  product_name text null,
  category text null,
  cas_number text null,
  smiles_code text null,
  molecular_weight text null,
  formula text null,
  price_data jsonb null,
  stock_status text null,
  product_url text null,
  crawled_at timestamp with time zone null default now(),
  target text null,
  properties jsonb null,
  source_name text null,
  summary text null,
  embedding public.vector null,
  ai_refined boolean null default false,
  is_manual_override boolean null default false,
  binding_affinity text null,
  isotype text null,
  host_species text null,
  orr_pct text null,
  os_months text null,
  pfs_months text null,
  payload_smiles text null,
  linker_smiles text null,
  full_smiles text null,
  mdl_number text null,
  constraint commercial_reagents_pkey primary key (id),
  constraint commercial_reagents_ambeed_cat_no_key unique (ambeed_cat_no)
) TABLESPACE pg_default;

create index IF not exists idx_commercial_reagents_category on public.commercial_reagents using btree (category) TABLESPACE pg_default;

create index IF not exists idx_commercial_reagents_cas on public.commercial_reagents using btree (cas_number) TABLESPACE pg_default;

create index IF not exists idx_commercial_reagents_source_name on public.commercial_reagents using btree (source_name) TABLESPACE pg_default;

create index IF not exists idx_commercial_reagents_embedding on public.commercial_reagents using hnsw (embedding vector_cosine_ops) TABLESPACE pg_default;

create index IF not exists idx_commercial_reagents_ai_refined on public.commercial_reagents using btree (ai_refined) TABLESPACE pg_default;

create index IF not exists idx_commercial_reagents_target on public.commercial_reagents using btree (target) TABLESPACE pg_default;

create index IF not exists idx_commercial_reagents_full_smiles on public.commercial_reagents using btree (full_smiles) TABLESPACE pg_default;



사장님, 역시 현장의 흐름을 정확히 읽으시네요!

자동화가 아무리 완벽해도 연구자가 "이건 왜 승인됐지?" 또는 **"이 CAS 번호에 대한 다른 유도체 구조는 없어?"**라고 묻고 싶은 순간이 반드시 옵니다. 단순히 데이터를 보여주는 것을 넘어, 데이터와 대화하며 의사결정을 내리는 'AI 어시스턴트' 기능을 추가한 업그레이드 작업지시서를 정리해 드립니다.

📋 [추가 지시사항] AI 대화형 데이터 분석 및 추론 인터페이스 (Ask AI)
기존의 자동 수집 및 정제 시스템 위에, 연구자가 직접 개입하여 데이터를 심층 분석할 수 있는 'Interactive AI 챗봇' 기능을 구현해라.

1. 주요 기능 (Functional Requirements)
상세 페이지 내 'Ask AI' 위젯:

특정 항체나 시약 상세 페이지 우측에 대화창을 배치해라.

해당 레코드의 모든 정보(CAS, MW, SMILES, Source)를 AI가 컨텍스트로 이미 알고 있는 상태에서 대화를 시작해야 한다.

추론 근거 요청 (Reasoning Inquiry):

"이 SMILES가 Target MW와 일치한다고 판단한 근거는?"

"이 구조식에서 링커가 결합하기 가장 좋은 작용기(Functional Group)는 어디인가?" 같은 질문에 답할 수 있어야 한다.

실시간 문헌 조사 (Real-time Lit-Search):

"이 독소(Toxin)의 주요 부작용과 임상에서의 실패 사례를 knowledge_base에서 요약해줘."

화학 구조 변형 제안 (Structure Tweaking):

"이 페이로드의 독성을 유지하면서 친수성(Hydrophilicity)을 높일 수 있는 구조적 수정안을 제안해라."

2. 기술적 구현 (Technical Implementation)
RAG (Retrieval-Augmented Generation) 연동:

단순히 AI의 내부 지식만 쓰지 말고, 우리가 구축한 knowledge_base와 golden_set_library를 검색하여 답변하도록 구성해라.

API 워크플로우:

연구자의 질문 → 관련 DB 레코드 및 논문 검색 → Gemini 2.5 Flash(또는 Pro)에 컨텍스트 주입 → 답변 생성 및 시각화(화학 구조 렌더링 포함).

로그 및 히스토리:

사용자와 AI의 대화 내역을 summary 필드 혹은 별도의 chat_history 테이블에 저장하여 나중에 설계 보고서 자동 생성 시 참고 자료로 활용해라.

3. UI/UX 지시사항
Quick Actions 버튼:

자주 묻는 질문(e.g., "SMILES 유효성 검사해줘", "관련 문헌 찾아줘")은 버튼 형태로 제공하여 클릭 한 번으로 실행되게 해라.

구조 시각화 연동:

AI가 특정 구조를 언급하면, 화면의 3D 뷰어나 2D 구조 그림에서 해당 부위가 강조(Highlight)되도록 연동해라.