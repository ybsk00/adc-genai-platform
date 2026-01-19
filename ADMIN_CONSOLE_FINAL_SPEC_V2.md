🚀 Admin Page Overhaul - Frontend Master Plan (Final V2)
Objective: 어드민페이지 개편 전체내역.md 및 세부내역.md에 기반하여 관리자 페이지의 UI를 전면 개편합니다. Mock 데이터를 제거하고, 실제 백엔드 API와 연동하여 데이터 수집-검수-발행 및 AI 튜닝의 전 과정을 통제하는 컨트롤 타워를 구축합니다.

1. 🔄 Data Operations (데이터 운영)
위치: /admin/data-operations (Tab Layout)

1.1 Data Sources (데이터 소스 제어)
UI Components:

Source Cards: ClinicalTrials.gov, PubMed, Perplexity 등 소스별 카드.

Status Badge: Last Run 시간, Success/Fail 상태, New Records 수 표시.

Actions:

Run Now: 즉시 크롤링 트리거.

Settings: 스케줄 설정 모달 (Cron Format UI).

API: POST /api/scheduler/sync/{source}, GET /api/scheduler/status.

1.2 Staging Area (검수소)
UI Components:

Inbox Table: status='draft'인 항목 리스트.

Comparison View (Split): 좌측 Raw Data vs 우측 Extracted JSON 비교 화면.

Actions:

✅ Approve: 승인 시 Golden Set으로 이동.

❌ Reject: 반려 (사유 입력 모달 -> 데이터 삭제).

✏️ Edit: AI 추출값 수정 후 승인.

API: GET /api/admin/drafts, POST .../approve, POST .../reject.

1.3 Golden Set Library (골든셋 관리)
UI Components:

Status Filter Tabs: All | Approved | Failed | Ongoing (탭으로 구분 필수).

Grid View: 약물명, 타겟, 링커, Outcome(성공/실패), IP Status(신호등) 컬럼.

Structure Viz: SMILES 코드를 smiles-drawer로 렌더링한 2D 이미지 컬럼.

Actions:

Manual Entry: + Add Data 버튼 -> 모달 폼. (SMILES 입력 시 하단에 실시간 구조 Preview 필수).

Bulk Upload: Excel/PDF 파일 Drag & Drop 업로드 (RAG 메타데이터 태깅).

API: GET /api/library/goldenset, POST /api/admin/goldenset.

1.4 Knowledge Base (지식 베이스)
UI Components:

List View: 뉴스/논문 리스트.

Evidence Card: Relevance Score, Source Tier (1~3), Key Facts (3줄 요약).

Reasoning Tooltip: Score에 마우스 오버 시 AI의 추천 사유("경쟁사 임상 실패 관련 중요 기사") 표시.

Actions:

Approve (Index): RAG 벡터 DB에 반영.

Trash: 삭제.

Bulk Action: ⚡ Auto-Approve Tier 1 (신뢰도 높은 소스 일괄 승인 버튼).

API: GET /api/knowledge-base, POST /api/knowledge-base/index.

2. 👤 User Operations (사용자 운영)
위치: /admin/user-operations

2.1 Simulation Logs (Admin View)
UI Components:

Monitoring Table: User, Project, Duration, Status, Error Summary.

Error Detail Modal: 실패 항목 클릭 시 전체 Stack Trace 표시.

Actions:

Retry: 관리자 권한으로 재실행 (User Credit 차감 없음).

Inspect: 결과 대시보드 열기.

API: GET /api/admin/simulations.

2.2 Credit Management & User List
UI Components:

User List: 검색, 플랜 정보, 잔여 크레딧.

Actions:

Grant/Revoke Credits: 크레딧 수동 조정.

Ban User: 악성 사용자 차단.

3. 🧠 AI Tuning (AI 튜닝)
위치: /admin/ai-tuning (3-Step Tabs)

3.1 Orchestrator (Main Router)
UI: 사용자 의도 파악 및 라우팅 로직을 위한 시스템 프롬프트 에디터.

3.2 Specialized Agents (Workers)
UI:

Agent Selector: Dropdown (Structure, Toxicology, Patent, Clinical).

Prompt Editor: 페르소나 및 CoT(Chain of Thought) 설정.

Test Sandbox: 프롬프트 수정 후 즉시 테스트 가능한 채팅창.

Version History: 이전 버전 복구 기능.

3.3 RAG Generator (Writer)
UI: 최종 답변 생성 스타일 및 인용구 포맷팅 지침 에디터.

API: GET /api/admin/prompts, PUT /api/admin/prompts/{agent_id}.

4. 🏠 Admin Overview (Dashboard Home)
위치: /admin/overview

System Status: API Server, Background Worker, DB 연결 상태 신호등.

Today's Metrics: 수집된 문서 수, 시뮬레이션 요청 수 KPI 카드.

Charts: 주간 사용자 증가 추이, 인기 타겟(Top Targets) 파이 차트.

Live Feed: 실시간 에러 및 중요 이벤트 로그 피드.

5. 🧪 My Projects (User View - New Page)
위치: /dashboard/my-projects (User-facing)

UI Components:

Project List: 날짜, 타겟, 상태 배지 (Processing - 애니메이션 / Completed / Failed).

Home Widgets: ADC Trend Feed (지식 베이스 기반 뉴스 요약), Quick Start.

Actions:

View Result: 상세 분석 페이지로 이동.

Download Report: PDF 리포트 다운로드.

Real-time: Supabase Realtime을 통한 상태 자동 업데이트.

⚠️ Backend Implementation Guidelines (Strict)
No New Microservices: 기존 backend/app 구조를 유지합니다.

Router Reuse: data_ops.py 등 새로운 파일을 무분별하게 생성하지 말고, 기존 admin.py 또는 data_processing.py에 로직을 추가(Append)하십시오.

Schema Updates: 데이터베이스 변경 시 DROP 하지 말고 반드시 **ALTER TABLE**을 사용하여 기존 데이터를 보존하십시오.
