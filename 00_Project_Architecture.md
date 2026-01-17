이 모든 파일의 기준이 되는 **00_Project_Architecture.md (전체 설계도)**를 작성해 드립니다. 이 문서를 보시면 제가 앞으로 작성할 나머지 9개 파일의 **디테일 수준(Depth)**을 가늠하실 수 있을 것입니다.

📄 00_Project_Architecture.md
Markdown

# 0. Project Overview & Architecture Design
**Project Name:** ADC-GenAI-Platform (Temporary Name)
**Version:** v1.0.0
**Date:** 2026-01-17
**Target Audience:** Bio-Researchers, CDMOs, VC Investment Reviewers

---

## 1. System Architecture Diagram (High Level)

[User Browser]
      │
      │ (HTTPS / JSON)
      ▼
[Frontend: React Single Page Application]
  │  • UI Library: Shadcn/UI + Tailwind CSS
  │  • State Management: TanStack Query (React Query)
  │  • Visualization: MolStar (3D), Chart.js (Radar)
  │
  ▼
[Backend: API Gateway (Python FastAPI)]
  │  • Auth: JWT / OAuth2 (Google, LinkedIn)
  │  • Payment: Lemon Squeezy Webhook Handler
  │
  ├──► [Database Layer]
  │       • PostgreSQL (User, Project, Payment Data)
  │       • pgvector (Vector Store for RAG)
  │       • Redis (Job Queue & Cache)
  │
  ▼
[AI Orchestration Layer (LangGraph)] ⚡ *The Brain*
  │
  ├──► **Controller Node (Orchestrator)**
  │       │
  │       ├──► Agent 1: Bio-Structure (BioNeMo/AlphaFold)
  │       ├──► Agent 2: Toxicology (RDKit + RAG)
  │       ├──► Agent 3: Patent (Search API)
  │       ├──► Agent 4: Competitor (Perplexity API)
  │       ├──► Agent 5: Clinical Design (LLM Reasoning)
  │       └──► Agent 6: Report Writer (LLM -> JSON)
  │
  ▼
[Output Generation]
  • PDF Generator (React-PDF Renderer)

---

## 2. Tech Stack Strategy (개발 스택 선정)

### 2.1 Frontend (Client Side)
* **Framework:** **React 18+ (Vite)** - *빠른 빌드와 렌더링 속도.*
* **Language:** **TypeScript** - *복잡한 바이오 데이터 타입 에러 방지.*
* **Styling:** **Tailwind CSS** - *빠른 UI 개발.*
* **UI Components:** **Shadcn/UI** - *깔끔하고 전문적인 Enterprise 룩앤필 제공.*
* **3D Viewer:** **Mol* (MolStar)** - *AlphaFold PDB 파일을 웹에서 가장 가볍게 돌리는 라이브러리.*
* **PDF Generation:** **React-PDF** - *화면 UI 그대로 PDF로 굽기 위해 사용.*

### 2.2 Backend (Server Side)
* **Framework:** **Python FastAPI** - *AI 라이브러리(LangChain, BioNeMo)와 호환성 최고.*
* **Task Queue:** **Celery + Redis** - *AI 분석이 3분 이상 걸리므로, 비동기(Async) 처리 필수.*
* **Orchestration:** **LangGraph** - *순차적 실행뿐만 아니라, 에이전트 간의 루프(Loop)와 조건부 분기 처리에 최적화.*

### 2.3 Database & Storage
* **Primary DB:** **Supabase (PostgreSQL)** - *관리 포인트 최소화, pgvector 플러그인 기본 지원.*
* **File Storage:** **AWS S3** (or Supabase Storage) - *생성된 PDB 파일과 PDF 리포트 저장.*

---

## 3. Directory Structure (폴더 구조 설계)

개발자가 프로젝트를 시작할 때 이 구조대로 폴더를 만듭니다.

```bash
root/
├── frontend/                 # React Application
│   ├── src/
│   │   ├── assets/           # Images, Fonts
│   │   ├── components/       # Reusable UI (Buttons, Inputs, Modals)
│   │   ├── pages/            # Dashboard, Builder, Admin
│   │   ├── hooks/            # Custom Hooks (useSimulation, useAuth)
│   │   ├── services/         # API Calls (Axios instances)
│   │   └── types/            # TypeScript Interfaces (Drug, Antibody)
│   └── ...
│
├── backend/                  # FastAPI Application
│   ├── app/
│   │   ├── api/              # API Endpoints (Routes)
│   │   ├── core/             # Config, Security (JWT)
│   │   ├── db/               # Database Models & Schemas
│   │   ├── services/         # Business Logic
│   │   │   ├── rag_service/  # Vector Search Logic
│   │   │   └── payments/     # Lemon Squeezy Integration
│   │   └── agents/           # 🤖 LangGraph Agents
│   │       ├── orchestrator.py
│   │       ├── structure_agent.py
│   │       ├── tox_agent.py
│   │       └── ...
│   └── worker/               # Celery Worker (Background Jobs)
│
├── data_pipeline/            # Data Collection Scripts
│   ├── crawlers/             # ClinicalTrials.gov Scraper
│   ├── parser/               # PDF to Chunk Parser
│   └── golden_set/           # Initial Excel Data
│
└── docs/                     # Design Documents (MD files)
4. Key Workflows (핵심 로직 흐름)
4.1 시뮬레이션 요청 흐름 (Synchronous vs Asynchronous)
사용자가 버튼을 클릭했을 때 브라우저가 멈추지 않도록 비동기(Polling) 방식을 사용합니다.

[Click] 사용자: "Start Simulation" 버튼 클릭.

[Request] 프론트엔드: 백엔드로 POST /api/simulation 요청 (항체 서열 포함).

[Response] 백엔드: 즉시 job_id: "adc_12345"를 반환하고 연결 끊음 (Non-blocking).

[Background] 백엔드 워커: AI 에이전트 6개 가동 시작.

[Polling] 프론트엔드: 5초마다 GET /api/simulation/{job_id}/status 요청.

Status: "Processing" (진행률 30%) → 로딩 바 업데이트.

Status: "Completed" → 결과 페이지로 이동.

5. Security & Compliance (보안 정책)
데이터 보안: 모든 항체 서열 데이터는 DB 저장 시 AES-256 암호화 처리 (고객 자산 보호).

API 보안: Rate Limiting 적용 (분당 60회 요청 제한) - API 과금 폭탄 방지.


---

### 💡 작성자의 코멘트 (다음 단계)

이 `00_Project_Architecture.md` 파일은 전체 숲을 보는 지도입니다.
이제 사용자님이 원하시는 **"버튼 클릭 후의 디테일"**은 다음 파일인 **`01_Frontend_Landing_Auth.md`**와 **`02_Frontend_User_Dashboard.md`**에서 본격적으로 나옵니다.

예를 들어, 다음 파일(`02`)에서는 아래와 같은 디테일이 들어갑니다.

> **[UI Component: Simulation Start Button]**
> * **Location:** 우측 하단 Floating Action Button.
> * **State:**
>     * `Default`: "Run Simulation" (Blue #007AFF)
>     * `Loading`: "Initializing Agents..." (Spinner 아이콘 회전)
>     * `Error`: "Failed (Retry)" (Red #FF3B30)
> * **Interaction:** 클릭 시 `useSimulationMutation` 훅 실행. `onSuccess` 시 Toast 메시지 "분석이 시작되었습니다" 노출.

**이 정도 깊이(Depth)로 `01`번 파일(프론트엔드 상세) 작성을 시작해도 될까요?**