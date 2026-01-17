Dashboard Layout Structure (레이아웃)
로그인 후 진입하는 모든 화면의 공통 레이아웃입니다.

코드 스니펫

graph TD
    A[Side Navigation Bar (LNB)] --> B[Top Header]
    B --> C[Main Content Area]
    A --> D[User Profile & Settings]
    C --> E[Dashboard Home]
    C --> F[ADC Builder (Simulation)]
    C --> G[My Projects (History)]
1.1 Side Navigation (Left Sidebar)
Width: 240px (Collapsible)

Menu Items:

🏠 Dashboard: 뉴스 피드, 크레딧 현황.

🧪 ADC Builder (New Project): [Core Feature] 시뮬레이션 실행.

📂 Library (Golden Set): 기존 데이터 검색.

📑 My Reports: 완료된 PDF 목록.

⚙️ Settings: API 키 관리, 팀원 초대.

1.2 Top Header
Breadcrumbs: 예: Home > ADC Builder > LIV-1 Project

Status Bar:

Credit: 🪙 450 Credits Available (클릭 시 충전 페이지 이동).

Server Status: 🟢 System Operational (Cloud Run 헬스체크 연동).

2. Page: Dashboard Home (메인 홈)
연구원이 출근해서 처음 보는 화면입니다.

Welcome Message: "Good Morning, Dr. [Name]. Ready to discover?"

Quick Actions (Cards):

[+] New Simulation (가장 크게 강조)

[📂] Browse Golden Set

Recent Activity (Table):

최근 돌렸던 시뮬레이션 5건의 상태 (Processing, Done, Failed).

ADC Trend Feed (RAG Widget):

"오늘의 ADC 뉴스: 화이자, 새로운 링커 기술 도입 발표..."

(백엔드 에이전트가 매일 아침 크롤링한 요약 정보를 보여줌).

3. Page: ADC Builder (핵심 기능 - 시뮬레이션 입력)
이 페이지가 서비스의 알파이자 오메가입니다. 3-Step Wizard 형태로 설계합니다.

Step 1. Target & Antibody (항체 입력)
Input Type:

Dropdown: 유명한 항체 선택 (Trastuzumab, Sacituzumab 등).

Manual Input: FASTA Sequence 직접 입력 (대형 텍스트 영역).

File Upload: .pdb 또는 .fasta 파일 업로드.

Validation:

입력 즉시 백엔드 API (POST /api/validate/sequence)를 호출하여 유효한 서열인지 체크.

Error: "Invalid Amino Acid character found at line 3."

Step 2. Payload & Linker (약물 접합)
UI: 좌측에는 옵션 선택, 우측에는 실시간 화학 구조(2D Structure) 미리보기 표시.

Payload Select:

Categories: Microtubule Inhibitors (MMAE), Topo1 Inhibitors (DXd), DNA Damagers.

Linker Select:

Options: Cleavable (Val-Cit), Non-cleavable (MCC), Custom SMILES Input.

DAR (Drug-to-Antibody Ratio):

Slider: 1 ~ 8 설정 (기본값 4).

Step 3. Configuration & Run (설정 및 실행)
Simulation Mode:

Fast Scan (1 Credit): 3D 구조만 빠르게 확인.

Deep Analysis (10 Credits): 독성, 특허, 경쟁사 분석 포함 (6-Agent 풀가동).

Job Name: 프로젝트 이름 입력 (예: LIV-1_MMAE_Test_01).

Action Button:

🚀 Run Simulation

클릭 시 로딩 애니메이션 → "Job Submitted" 토스트 메시지 → Result Page로 자동 이동.

4. Page: Result Viewer (결과 화면)
백엔드(Cloud Run)가 열심히 계산하는 동안, 그리고 계산이 끝난 후 보여주는 화면입니다.

4.1 Loading State (Progress View)
Polling Logic: 프론트엔드는 5초마다 GET /api/simulation/{id}/status를 호출합니다.

Visual Stepper: 6명의 에이전트 작업 현황을 시각적으로 보여줍니다.

[✓] Structure Agent ... Done

[↻] Toxicity Agent ... Running (45%)

[ ] Patent Agent ... Pending

4.2 Completed View (Report Dashboard)
PDF를 다운로드하기 전, 웹에서 핵심 결과를 먼저 보여줍니다.

Summary Card:

Grade: B+ (Color Coded)

Verdict: "Go for In-vitro"

3D Viewer (MolStar Integration):

중앙에 ADC 3D 모델 렌더링.

마우스 조작(회전, 확대) 가능.

기능: Export PDB 버튼 제공.

Charts:

Radar Chart: [효능, 독성, 물성, 특허, 생산성] 5각 편대 분석.

Action:

📥 Download Full PDF Report

Method,Endpoint,Description,Request Body
POST,/api/jobs,시뮬레이션 시작 요청,"{ ""antibody_seq"": ""..."", ""payload_id"": ""mmae"", ""mode"": ""deep"" }"
GET,/api/jobs/{id},진행 상태 및 결과 조회,-
GET,/api/library/goldenset,골든셋 목록 검색,?target=LIV-1&page=1