03. Frontend Design: Admin Console
Document ID: FE-03 Role: Business Operation, Data Management, AI Tuning Access Control: Role == 'admin' Only (Strict RLS Enforcement) Tech Stack: React, Shadcn/UI (Data Table), Recharts

1. Admin Layout & Security (보안 및 구조)
1.1 Route Protection (접근 제어)
URL: /admin/*

Logic:

페이지 진입 시 useAuth() 훅을 통해 사용자의 role을 확인.

admin이 아니면 즉시 메인 홈(/)으로 리다이렉트 (403 Forbidden).

1.2 Layout Structure
Color Theme: 일반 사용자용(파란색/흰색)과 구분하기 위해 **다크 모드(Slate-900)**를 기본으로 적용하여 "관리자 모드임"을 시각적으로 강조.

Admin Sidebar:

📊 Overview: 핵심 KPI 대시보드.

👥 User Ops: 회원 관리 및 크레딧 지급.

🧬 Data Ops: 골든셋(Golden Set) 데이터 수정.

🤖 AI Tuning: 에이전트 프롬프트 관리.

⚙️ System: 서버 상태 및 결제 로그.

2. Page: Admin Overview (현황판)
사장님이 출근해서 커피 마시며 보는 화면입니다.

KPI Cards (Top Row):

MRR (Monthly Recurring Revenue): $12,450 (이번 달 예상 매출).

Total Users: 1,204 (신규 +15).

Active Simulations: 45 Jobs running now.

Error Rate: 0.5% (지난 24시간 기준).

Charts:

User Growth: 주간 가입자 추이 (Line Chart).

Top Targets: 사용자들이 가장 많이 분석한 타겟 (예: LIV-1 40%, HER2 30%) - 시장 트렌드 파악용.

Real-time Feed:

"User kim@bionet.com purchased Pro Plan ($499)."

"Simulation #Job-992 failed (Timeout)."

3. Page: User Operations (회원 관리)
고객 불만 처리(CS)와 영업 지원을 위한 페이지입니다.

3.1 User List Table
Columns: Name, Email, Organization, Plan, Credits, Last Login.

Search: 이메일 또는 소속으로 검색.

Quick Actions (Row Menu):

View Details: 해당 유저가 돌린 시뮬레이션 기록 열람.

Ban User: 악성 유저 정지.

💰 Grant Credits (핵심): 수동으로 크레딧 넣어주기.

3.2 Feature: Manual Credit Injection (크레딧 지급 모달)
Scenario: 영업차 만난 'Syngene' 팀장에게 체험판 포인트를 주고 싶음.

Modal UI:

Target User: lee@syngene.com

Amount: [ 500 ] Credits

Reason: 영업용 체험판 지급 (Memo)

[Send] Button: 클릭 시 즉시 반영되고, 유저에게 "관리자가 500 크레딧을 선물했습니다" 이메일 발송.

4. Page: Data Ops (Golden Set Manager)
개발자에게 부탁하지 않고, 대표님이 직접 엑셀처럼 데이터를 수정하는 곳입니다.

4.1 Golden Set Editor (Data Grid)
UI: 엑셀과 똑같은 UI (Ag-Grid 또는 TanStack Table).

Function:

Cell Edit: '엔허투'의 임상 단계를 'Phase 2' → 'Approved'로 더블클릭해서 수정.

Add Row: 새로운 약물 데이터 추가.

Import CSV: 엑셀 파일 통째로 업로드하여 DB 덮어쓰기.

RAG Trigger:

[🔄 Re-index Vector DB] Button: 데이터를 수정했으면 이 버튼을 눌러야 AI(RAG)가 바뀐 내용을 공부합니다. (필수 기능)

5. Page: AI Tuning (Prompt Manager)
AI가 헛소리를 하거나, 리포트 말투가 마음에 안 들 때 코드를 건드리지 않고 프롬프트만 수정합니다.

Agent Selector: 탭(Tab)으로 구분 (Structure Agent, Toxicology Agent, Report Writer).

Prompt Editor (Text Area):

Current System Prompt:

"You are a toxicology expert. Analyze the risks..."

Modification:

"... Analyze the risks and explicitly mention 'Ocular Toxicity' if the linker is unstable." (강조 구문 추가)

Test Sandbox:

프롬프트 수정 후 바로 테스트해 볼 수 있는 채팅창 제공.

"LIV-1 테스트해 봐." → 결과 확인 → [Save & Deploy] 버튼 클릭.

6. API Requirements (Admin Only)
일반 API와 달리 강력한 권한을 가진 엔드포인트입니다.

Method,Endpoint,Description,Request Body
GET,/api/admin/stats,대시보드 KPI 조회,-
POST,/api/admin/credits/grant,특정 유저 크레딧 지급,"{ ""user_id"": ""..."", ""amount"": 500, ""reason"": ""Sales"" }"
POST,/api/admin/goldenset/sync,벡터 DB(RAG) 재동기화,-
PUT,/api/admin/prompts/{agent_id},에이전트 프롬프트 수정,"{ ""new_prompt"": ""..."" }"

7. Development Checklist
[ ] Admin Layout (Sidebar, Dark Mode) 구현.

[ ] User Table 및 Credit Grant Modal 구현.

[ ] Prompt Manager UI 구현 (DB의 prompts 테이블과 연동).

[ ] RLS(Row Level Security) 확인 (일반 유저가 API 호출 시 차단되는지 테스트).