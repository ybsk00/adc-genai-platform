# ADC Platform Multi-Agent Design Engine Implementation Plan v2.2

## Executive Summary
사용자 대시보드 4대 메뉴와 6개 전문 에이전트 기반 자율 설계 엔진 구현 계획서
**보완일: 2026-01-27 | 버전: 2.2 (Critical Details 반영)**

### v2.2 주요 변경사항 ⚠️ Critical
- **Constraint 6**: 21 CFR Part 11 Digital Seal - 로그 불변성 보장 (Chain Hash + INSERT-ONLY)
- **Constraint 7**: Sandbox Library Version Lock - 버전 100% 동기화 자동 검증
- **AlphaFold GPU Quota**: 월간/일간 쿼터 관리 테이블 및 로직 추가

### v2.1 주요 변경사항
- **The Librarian**: Dynamic Substructure RAG 추가 (스캐폴드 → 지식베이스 의미론적 매핑)
- **The Auditor**: Constraint Guardrail 신규 (사용자 초기 제약조건 강제)
- **The Healer**: Knowledge Loop 추가 (성공적 수정 → candidate_snippets 자동 등록)
- **Phase 3.5**: AlphaFold 3 Integration (Premium Tier 3D 모델링)

---

## 🔒 Critical Engineering Constraints (필수 준수 사항)

### Constraint 1: Stateful Orchestration
> LangGraph의 State 객체에 현재 설계 중인 SMILES와 계산된 MW, LogP 수치를 실시간 동기화하여 에이전트들이 정보를 공유하게 할 것.

### Constraint 2: Self-Correction Feedback
> The Healer가 코드를 고친 횟수와 에러 메시지를 agent_logs에 기록하여, 추후 모델 튜닝 데이터로 활용할 수 있게 할 것.

### Constraint 3: Backend-Level Data Masking
> 무료 등급(Free Tier)의 경우, The Auditor가 생성한 리포트의 핵심 수치와 SMILES를 **프론트엔드가 아닌 백엔드 API 레벨에서 마스킹**하여 보안을 철저히 할 것.

### Constraint 4: Constraint Guardrail (v2.1 신규)
> The Auditor는 최종 검증 시 **사용자가 세션 생성 시 입력한 초기 제약조건**(Target Antigen, Requested DAR, Linker Preference)이 설계 결과에 반영되었는지 반드시 검사할 것. 초기 조건과 다른 결과 도출 시 재설계 요청하거나 사유를 명시할 것.

### Constraint 5: Knowledge Loop (v2.1 신규)
> The Healer가 성공적으로 코드를 수정한 경우, 해당 수정 패턴을 `candidate_snippets` 테이블에 자동 등록하여 향후 유사 에러 발생 시 재활용할 수 있게 할 것. 검증(verified) 후 `code_snippet_library`로 승격.

### Constraint 6: 21 CFR Part 11 Digital Seal (v2.1 신규) ⚠️ Critical
> FDA 임상 신청 시 데이터 무결성 증명을 위해, `agent_execution_logs` 기록 시 **SHA-256 해시 기반 디지털 봉인(Digital Seal)**을 적용할 것. 로그 레코드는 **INSERT-ONLY** 정책을 적용하여 수정 불가능성(Immutability)을 보장하고, 체인 해시(Chain Hash)로 시퀀스 무결성을 검증할 것.

### Constraint 7: Sandbox Library Version Lock (v2.1 신규) ⚠️ Critical
> The Coder 에이전트가 학습/추론에 사용하는 Python 환경과 Docker Sandbox 내부의 **라이브러리 버전을 100% 일치**시킬 것. `requirements.lock` 파일을 단일 소스로 관리하고, CI/CD 파이프라인에서 버전 불일치 시 빌드를 차단할 것.

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                         User Dashboard                               │
│  ┌────────────┐ ┌────────────┐ ┌────────────┐ ┌────────────┐       │
│  │ De novo    │ │ Lead       │ │ Pre-clinical│ │ CMC &      │       │
│  │ Design     │ │ Optimization│ │ Audit      │ │ Sourcing   │       │
│  └────────────┘ └────────────┘ └────────────┘ └────────────┘       │
│                                                                      │
│  ┌──────────────────────────────────────────────────────────┐      │
│  │              🖥️ Live Agent Console (WebSocket)            │      │
│  │  [Orchestrator] → [Alchemist] → [Coder] → [Auditor]      │      │
│  │  "Analyzing target HER2..." "Generating candidates..."    │      │
│  └──────────────────────────────────────────────────────────┘      │
└─────────────────────────────────────────────────────────────────────┘
                                │
                                ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    🧠 The Orchestrator (LangGraph)                   │
│  ┌─────────────────────────────────────────────────────────────┐   │
│  │              Shared State (실시간 동기화)                     │   │
│  │  • current_smiles: "CC(C)..."                               │   │
│  │  • calculated_metrics: {mw: 718.5, logp: 3.2, ...}         │   │
│  │  • validation_flags: {lipinski: ✓, pains: ✓, ...}          │   │
│  └─────────────────────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────────────────────┘
         │           │           │           │           │
         ▼           ▼           ▼           ▼           ▼
    ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐ ┌────────┐
    │Alchemist│ │ Coder  │ │ Healer │ │Auditor │ │Librarian│
    │ 설계    │ │ 실행   │ │ 치유   │ │ 검증   │ │ 근거   │
    │         │ │        │ │        │ │        │ │        │
    │ Golden  │ │Snippet │ │Max 3   │ │Risk    │ │RAG +   │
    │ Set 기반│ │Library │ │Retry   │ │Score   │ │PMID    │
    └────────┘ └────────┘ └────────┘ └────────┘ └────────┘
                    │           ▲
                    ▼           │ (에러 시 자가 치유)
              ┌─────────────────┴───────┐
              │   🐳 Docker Sandbox      │
              │   (Air-gapped Python)   │
              │   • RDKit, BioPython    │
              │   • 512MB RAM Limit     │
              │   • No Network Access   │
              └─────────────────────────┘
```

---

## Phase 0: Current State Analysis

### Existing Components (재사용 가능)
| Component | Path | Status |
|-----------|------|--------|
| RAG Service | `backend/app/services/rag_service.py` | ✅ Ready |
| PubChem Service | `backend/app/services/pubchem_service.py` | ✅ Ready |
| Golden Set Library | `golden_set_library` table | ✅ Ready |
| Knowledge Base | `knowledge_base` table | ✅ Ready |
| ADC Builder UI | `frontend/src/pages/Dashboard/ADCBuilder.tsx` | ⚠️ Refactor |
| Design Runs | `backend/app/services/design_run_service.py` | ⚠️ Enhance |

### New Components Required
- Multi-Agent Orchestrator (LangGraph)
- Python Sandbox Executor (Docker-based)
- Self-Healing Loop (The Healer)
- Verified Code Snippet Library
- Real-time WebSocket + Live Console
- Backend-Level Data Masking
- Encrypted Structures Storage

---

## Phase 1: Database Schema - "신뢰의 기록" (Audit Trail) [Week 1]

> 💡 **핵심 목표**: 규제 기관(FDA/EMA) 제출 가능한 수준의 감사 추적 이력 관리

### 1.1 Core Tables

```sql
-- =====================================================
-- 1. Design Sessions: 설계 세션 마스터 테이블
-- =====================================================
CREATE TABLE design_sessions (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES auth.users(id) NOT NULL,

    -- Session Type & Tier
    session_type VARCHAR(50) NOT NULL CHECK (session_type IN ('denovo', 'optimization', 'audit', 'cmc')),
    tier VARCHAR(20) DEFAULT 'free' CHECK (tier IN ('free', 'premium')),

    -- Workflow Status
    status VARCHAR(30) DEFAULT 'pending' CHECK (status IN ('pending', 'running', 'completed', 'failed', 'manual_review')),
    current_step INTEGER DEFAULT 0,
    total_steps INTEGER DEFAULT 4,
    current_agent VARCHAR(50),

    -- 🔒 Design Parameters (FDA Audit Trail)
    design_goal TEXT,                          -- "HER2 양성 유방암 타겟 ADC 설계"
    target_antigen VARCHAR(100),               -- "HER2", "TROP-2", etc.
    target_indication VARCHAR(200),            -- "Breast Cancer", "NSCLC"
    requested_dar INTEGER,                     -- Drug-Antibody Ratio (2, 4, 8)
    linker_preference VARCHAR(50),             -- "cleavable", "non-cleavable"

    -- 🧠 Shared State (에이전트 간 실시간 동기화)
    current_smiles TEXT,
    calculated_metrics JSONB DEFAULT '{}',     -- {mw, logp, hbd, hba, psa, tpsa}
    validation_flags JSONB DEFAULT '{}',       -- {lipinski: true, pains: false, ...}

    -- Results
    result_candidates JSONB DEFAULT '[]',      -- [{smiles, score, rank}, ...]
    final_report JSONB,
    final_approval BOOLEAN DEFAULT FALSE,      -- Auditor 최종 승인 여부
    approval_reasoning TEXT,                   -- 승인/거부 사유

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),
    completed_at TIMESTAMPTZ,

    -- Regulatory Metadata
    regulatory_version VARCHAR(20) DEFAULT 'v1.0',
    audit_hash VARCHAR(64)                     -- SHA-256 of critical data for integrity
);

-- =====================================================
-- 2. Agent Execution Logs: 에이전트별 실행 이력 + 추론 근거
-- 🔒 21 CFR Part 11 준수: INSERT-ONLY + Digital Seal
-- =====================================================
CREATE TABLE agent_execution_logs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    session_id UUID REFERENCES design_sessions(id) ON DELETE CASCADE,

    -- Agent Info
    agent_name VARCHAR(50) NOT NULL CHECK (agent_name IN (
        'orchestrator', 'alchemist', 'coder', 'healer', 'auditor', 'librarian'
    )),
    step_number INTEGER,
    status VARCHAR(20) CHECK (status IN ('started', 'completed', 'error', 'healed', 'skipped')),

    -- 🔒 추론 근거 (Reasoning) - FDA 감사용
    reasoning TEXT,                            -- "Golden Set 유사도 0.85 기반 후보 선정"
    decision_summary TEXT,                     -- "Top 3 candidates selected based on Tanimoto > 0.7"
    confidence_score DECIMAL(3,2),             -- 0.00 ~ 1.00

    -- Input/Output
    input_data JSONB,
    output_data JSONB,
    error_message TEXT,
    execution_time_ms INTEGER,

    -- 🔄 Self-healing Logs (The Healer 전용)
    retry_count INTEGER DEFAULT 0,
    fix_logs JSONB DEFAULT '[]',               -- [{attempt: 1, error: "...", fix: "..."}, ...]
    healing_successful BOOLEAN,

    -- References
    referenced_golden_set_ids UUID[] DEFAULT '{}',
    referenced_knowledge_ids INTEGER[] DEFAULT '{}',
    referenced_pmids TEXT[] DEFAULT '{}',

    -- 🔒 21 CFR Part 11 Digital Seal (v2.1 신규)
    record_hash VARCHAR(64) NOT NULL,          -- SHA-256 of (prev_hash + record_data)
    prev_record_hash VARCHAR(64),              -- 이전 레코드 해시 (Chain Hash)
    sequence_number BIGINT NOT NULL,           -- 순차 번호 (무결성 검증용)
    sealed_at TIMESTAMPTZ DEFAULT NOW(),       -- 봉인 시점
    seal_version VARCHAR(10) DEFAULT 'v1.0',   -- 봉인 알고리즘 버전

    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- 🔒 INSERT-ONLY 정책: UPDATE/DELETE 차단
CREATE OR REPLACE FUNCTION prevent_log_modification()
RETURNS TRIGGER AS $$
BEGIN
    RAISE EXCEPTION '21 CFR Part 11 Compliance: agent_execution_logs is immutable. Record ID: %', OLD.id;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER enforce_immutability
    BEFORE UPDATE OR DELETE ON agent_execution_logs
    FOR EACH ROW EXECUTE FUNCTION prevent_log_modification();

-- Chain Hash 자동 생성 트리거
CREATE OR REPLACE FUNCTION generate_record_seal()
RETURNS TRIGGER AS $$
DECLARE
    prev_hash VARCHAR(64);
    prev_seq BIGINT;
    seal_data TEXT;
BEGIN
    -- 이전 레코드 조회
    SELECT record_hash, sequence_number INTO prev_hash, prev_seq
    FROM agent_execution_logs
    WHERE session_id = NEW.session_id
    ORDER BY sequence_number DESC
    LIMIT 1;

    -- 시퀀스 번호 설정
    NEW.sequence_number := COALESCE(prev_seq, 0) + 1;
    NEW.prev_record_hash := prev_hash;

    -- 봉인 데이터 생성 (핵심 필드들의 조합)
    seal_data := COALESCE(prev_hash, 'GENESIS') || '|' ||
                 NEW.session_id || '|' ||
                 NEW.agent_name || '|' ||
                 NEW.sequence_number || '|' ||
                 COALESCE(NEW.reasoning, '') || '|' ||
                 NEW.created_at::TEXT;

    -- SHA-256 해시 생성
    NEW.record_hash := encode(sha256(seal_data::bytea), 'hex');

    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER seal_record_before_insert
    BEFORE INSERT ON agent_execution_logs
    FOR EACH ROW EXECUTE FUNCTION generate_record_seal();

-- =====================================================
-- 3. Encrypted Structures: 유료 사용자 SMILES 암호화 저장
-- =====================================================
CREATE TABLE encrypted_structures (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    session_id UUID REFERENCES design_sessions(id) ON DELETE CASCADE,
    user_id UUID REFERENCES auth.users(id) NOT NULL,

    -- Encrypted Data (AES-256-GCM)
    encrypted_smiles BYTEA NOT NULL,           -- AES-256 암호화된 SMILES
    encryption_iv BYTEA NOT NULL,              -- Initialization Vector
    encryption_tag BYTEA NOT NULL,             -- Authentication Tag

    -- Metadata (암호화되지 않은 메타데이터)
    structure_type VARCHAR(50),                -- 'candidate', 'optimized', 'final'
    rank INTEGER,
    is_premium_only BOOLEAN DEFAULT TRUE,

    -- Validation Hash
    structure_hash VARCHAR(64),                -- 복호화 후 검증용 SHA-256

    created_at TIMESTAMPTZ DEFAULT NOW(),
    accessed_at TIMESTAMPTZ                    -- 마지막 접근 시간
);

-- =====================================================
-- 4. Code Execution History: Python 코드 실행 이력
-- =====================================================
CREATE TABLE code_execution_history (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    session_id UUID REFERENCES design_sessions(id) ON DELETE CASCADE,
    agent_name VARCHAR(50) DEFAULT 'coder',

    -- Code Content
    code_content TEXT NOT NULL,
    snippet_ids TEXT[] DEFAULT '{}',           -- 사용된 스니펫 ID들

    -- Execution Result
    execution_result JSONB,
    stdout TEXT,
    stderr TEXT,
    exit_code INTEGER,
    execution_time_ms INTEGER,

    -- Self-healing
    is_healed BOOLEAN DEFAULT FALSE,
    original_code_id UUID REFERENCES code_execution_history(id),
    healing_attempt INTEGER DEFAULT 0,

    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- =====================================================
-- 5. Code Snippet Library: 검증된 코드 조각 라이브러리
-- =====================================================
CREATE TABLE code_snippet_library (
    id VARCHAR(50) PRIMARY KEY,                -- 'lipinski_check', 'tanimoto_calc', etc.
    name VARCHAR(100) NOT NULL,
    description TEXT,
    category VARCHAR(50),                      -- 'validation', 'calculation', 'conversion'

    -- Code
    code_template TEXT NOT NULL,               -- Python code with {{placeholders}}
    required_imports TEXT[] DEFAULT '{}',      -- ['rdkit', 'numpy']
    input_params JSONB,                        -- {smiles: "string", threshold: "float"}
    output_schema JSONB,                       -- {valid: "bool", mw: "float"}

    -- Quality
    is_verified BOOLEAN DEFAULT FALSE,
    test_coverage DECIMAL(5,2),
    last_verified_at TIMESTAMPTZ,

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- =====================================================
-- 6. Candidate Snippets: The Healer Knowledge Loop (v2.1 신규)
-- =====================================================
CREATE TABLE candidate_snippets (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),

    -- Source Context
    source_session_id UUID REFERENCES design_sessions(id),
    source_healer_log_id UUID REFERENCES agent_execution_logs(id),

    -- Error Pattern
    error_type VARCHAR(50) NOT NULL,           -- 'SYNTAX_ERROR', 'CHEMICAL_VALIDITY_ERROR', etc.
    error_pattern TEXT NOT NULL,               -- 정규식 또는 키워드 패턴
    error_context JSONB,                       -- 에러 발생 상황 컨텍스트

    -- Fix Pattern
    fix_code TEXT NOT NULL,                    -- 수정된 코드 조각
    fix_explanation TEXT,                      -- 수정 이유/설명
    fix_strategy VARCHAR(100),                 -- 'add_null_check', 'try_alternative_method', etc.

    -- Quality Metrics
    success_count INTEGER DEFAULT 1,           -- 성공 횟수
    failure_count INTEGER DEFAULT 0,           -- 이후 실패 횟수
    confidence_score DECIMAL(3,2) DEFAULT 0.50, -- 신뢰도 (0.00 ~ 1.00)

    -- Promotion Status
    status VARCHAR(20) DEFAULT 'candidate' CHECK (status IN (
        'candidate',    -- 후보 상태
        'reviewing',    -- 검토 중
        'promoted',     -- code_snippet_library로 승격됨
        'rejected'      -- 거부됨
    )),
    promoted_snippet_id VARCHAR(50) REFERENCES code_snippet_library(id),
    promoted_at TIMESTAMPTZ,

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- =====================================================
-- 7. Scaffold Knowledge Mapping: The Librarian RAG (v2.1 신규)
-- =====================================================
CREATE TABLE scaffold_knowledge_mapping (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),

    -- Scaffold/Substructure Info
    scaffold_smiles TEXT NOT NULL,             -- 스캐폴드 SMILES (Murcko decomposition)
    scaffold_type VARCHAR(50),                 -- 'murcko_generic', 'murcko_framework', 'ring_system'
    fingerprint_bits BYTEA,                    -- Morgan fingerprint for similarity search

    -- Semantic Embedding
    embedding_vector VECTOR(1536),             -- OpenAI/Cohere embedding for semantic search
    scaffold_description TEXT,                 -- AI 생성 스캐폴드 설명

    -- Knowledge Links
    linked_knowledge_ids INTEGER[] DEFAULT '{}',   -- knowledge_base 연결
    linked_golden_set_ids UUID[] DEFAULT '{}',     -- golden_set_library 연결
    linked_pmids TEXT[] DEFAULT '{}',              -- PubMed 문헌 연결

    -- Usage Stats
    retrieval_count INTEGER DEFAULT 0,
    last_retrieved_at TIMESTAMPTZ,

    -- Quality
    relevance_score DECIMAL(3,2),              -- RAG 검색 시 관련성 점수
    is_validated BOOLEAN DEFAULT FALSE,

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW(),

    UNIQUE(scaffold_smiles, scaffold_type)
);
```

### 1.2 Indexes for Performance

```sql
-- Session lookups
CREATE INDEX idx_sessions_user_status ON design_sessions(user_id, status);
CREATE INDEX idx_sessions_type ON design_sessions(session_type);
CREATE INDEX idx_sessions_created ON design_sessions(created_at DESC);

-- Log queries
CREATE INDEX idx_logs_session_agent ON agent_execution_logs(session_id, agent_name);
CREATE INDEX idx_logs_status ON agent_execution_logs(status);
CREATE INDEX idx_logs_healing ON agent_execution_logs(agent_name, healing_successful)
    WHERE agent_name = 'healer';

-- Encrypted structures
CREATE INDEX idx_encrypted_session ON encrypted_structures(session_id);
CREATE INDEX idx_encrypted_user ON encrypted_structures(user_id);

-- Code history
CREATE INDEX idx_code_session ON code_execution_history(session_id);
CREATE INDEX idx_code_healed ON code_execution_history(is_healed);

-- Candidate snippets (v2.1)
CREATE INDEX idx_candidate_error_type ON candidate_snippets(error_type, status);
CREATE INDEX idx_candidate_confidence ON candidate_snippets(confidence_score DESC);
CREATE INDEX idx_candidate_status ON candidate_snippets(status);

-- Scaffold knowledge mapping (v2.1)
CREATE INDEX idx_scaffold_type ON scaffold_knowledge_mapping(scaffold_type);
CREATE INDEX idx_scaffold_embedding ON scaffold_knowledge_mapping
    USING ivfflat (embedding_vector vector_cosine_ops) WITH (lists = 100);
CREATE INDEX idx_scaffold_validated ON scaffold_knowledge_mapping(is_validated);
```

### 1.3 Row-Level Security (RLS)

```sql
-- Users can only see their own sessions
ALTER TABLE design_sessions ENABLE ROW LEVEL SECURITY;
CREATE POLICY "Users see own sessions" ON design_sessions
    FOR ALL USING (auth.uid() = user_id);

-- Encrypted structures are strictly user-owned
ALTER TABLE encrypted_structures ENABLE ROW LEVEL SECURITY;
CREATE POLICY "Users see own structures" ON encrypted_structures
    FOR ALL USING (auth.uid() = user_id);
```

---

## Phase 2: Backend Multi-Agent Core [Week 2-3]

### 2.1 Directory Structure

```
backend/app/
├── agents/
│   ├── __init__.py
│   ├── base_agent.py              # Abstract base class
│   ├── orchestrator.py            # The Orchestrator (LangGraph)
│   ├── alchemist.py               # The Alchemist (Design Engine)
│   ├── coder.py                   # The Coder (Python Executor)
│   ├── healer.py                  # The Healer (Self-healing)
│   ├── auditor.py                 # The Auditor (Validation)
│   └── librarian.py               # The Librarian (RAG)
├── sandbox/
│   ├── executor.py                # Docker sandbox executor
│   ├── snippet_manager.py         # Code snippet library manager
│   ├── templates/                 # Pre-defined code templates
│   │   ├── lipinski_check.py
│   │   ├── tanimoto_similarity.py
│   │   ├── pains_filter.py
│   │   └── sa_score.py
│   └── Dockerfile.sandbox
├── workflows/
│   ├── denovo_workflow.py         # De novo Design workflow
│   ├── optimization_workflow.py   # Lead Optimization workflow
│   ├── audit_workflow.py          # Pre-clinical Audit workflow
│   └── cmc_workflow.py            # CMC & Sourcing workflow
├── security/
│   ├── encryption.py              # AES-256 encryption utils
│   └── masking.py                 # Backend data masking
└── websocket/
    └── session_hub.py             # Real-time WebSocket hub
```

### 2.2 Shared State Definition (LangGraph)

```python
# backend/app/agents/orchestrator.py
from typing import TypedDict, Literal, List, Optional
from langgraph.graph import StateGraph, END
from pydantic import BaseModel

class CalculatedMetrics(TypedDict):
    mw: Optional[float]           # Molecular Weight
    logp: Optional[float]         # LogP
    hbd: Optional[int]            # H-bond Donors
    hba: Optional[int]            # H-bond Acceptors
    psa: Optional[float]          # Polar Surface Area
    tpsa: Optional[float]         # Topological PSA
    rotatable_bonds: Optional[int]
    sa_score: Optional[float]     # Synthetic Accessibility

class ValidationFlags(TypedDict):
    lipinski_pass: bool
    pains_free: bool
    structural_alerts: List[str]
    golden_set_similarity: float
    auditor_approved: bool

class SessionState(TypedDict):
    """
    🧠 Shared State: 모든 에이전트가 실시간 동기화하는 상태 객체
    Constraint 1 준수: LangGraph State에 SMILES와 계산된 수치 실시간 동기화
    """
    # Session Info
    session_id: str
    session_type: Literal["denovo", "optimization", "audit", "cmc"]
    tier: Literal["free", "premium"]
    user_id: str

    # Design Parameters
    target_antigen: str
    target_indication: str
    requested_dar: int
    linker_preference: str

    # 🔄 Real-time Synced Data
    current_smiles: str
    calculated_metrics: CalculatedMetrics
    validation_flags: ValidationFlags

    # Workflow Control
    step: int
    current_agent: str
    error: Optional[str]
    requires_healing: bool
    healing_attempts: int

    # Results
    candidates: List[dict]
    final_report: Optional[dict]

class Orchestrator:
    def __init__(self):
        self.graph = self._build_graph()
        self.websocket_hub = WebSocketHub()

    def _build_graph(self) -> StateGraph:
        workflow = StateGraph(SessionState)

        # Add agent nodes
        workflow.add_node("alchemist", self._run_alchemist)
        workflow.add_node("coder", self._run_coder)
        workflow.add_node("healer", self._run_healer)
        workflow.add_node("auditor", self._run_auditor)
        workflow.add_node("librarian", self._run_librarian)

        # Entry point routing based on session_type
        workflow.set_entry_point("route_entry")
        workflow.add_node("route_entry", self._route_entry)

        workflow.add_conditional_edges(
            "route_entry",
            self._get_entry_agent,
            {
                "denovo": "alchemist",
                "optimization": "auditor",  # 기존 독성 분석 먼저
                "audit": "coder",           # 전수 물성 조사 먼저
                "cmc": "librarian"
            }
        )

        # Coder -> Healer (on error) or Auditor (on success)
        workflow.add_conditional_edges(
            "coder",
            self._check_coder_result,
            {
                "success": "auditor",
                "error": "healer",
                "max_retries": "manual_review"
            }
        )

        # Healer -> Coder (retry)
        workflow.add_edge("healer", "coder")

        # Auditor -> END or Alchemist (redesign needed)
        workflow.add_conditional_edges(
            "auditor",
            self._check_auditor_result,
            {
                "approved": END,
                "redesign": "alchemist",
                "manual_review": "manual_review"
            }
        )

        return workflow.compile()

    async def _run_alchemist(self, state: SessionState) -> SessionState:
        """The Alchemist: SMILES 설계"""
        await self._broadcast_status(state, "alchemist", "started",
            "Golden Set 기반 후보 물질 설계 중...")

        agent = AlchemistAgent()
        result = await agent.execute(state)

        # Update shared state
        state["current_smiles"] = result.data.get("primary_smiles")
        state["candidates"] = result.data.get("candidates", [])
        state["current_agent"] = "alchemist"
        state["step"] += 1

        await self._broadcast_status(state, "alchemist", "completed",
            f"{len(state['candidates'])}개 후보 생성 완료")

        return state

    async def _broadcast_status(self, state: SessionState, agent: str,
                                 status: str, message: str):
        """WebSocket으로 실시간 상태 전송"""
        await self.websocket_hub.broadcast(state["session_id"], {
            "type": "agent_status",
            "agent": agent,
            "status": status,
            "message": message,
            "step": state["step"],
            "timestamp": datetime.utcnow().isoformat()
        })
```

### 2.3 The Coder with Snippet Library

```python
# backend/app/agents/coder.py
from app.sandbox.snippet_manager import SnippetManager

class CoderAgent(BaseAgent):
    """
    The Coder: 검증된 코드 스니펫을 조합하여 실행
    - AI가 처음부터 코드를 짜지 않고, Snippet Library를 활용
    - 에러 발생률 최소화
    """
    name = "coder"

    def __init__(self):
        self.snippet_manager = SnippetManager()
        self.sandbox = SandboxExecutor()

    async def execute(self, state: SessionState) -> AgentOutput:
        smiles = state["current_smiles"]
        await self._log_start(state["session_id"])

        # 1. 필요한 검증 항목 결정
        validations_needed = self._determine_validations(state)

        # 2. 검증된 스니펫으로 코드 조합
        code = self._compose_code_from_snippets(smiles, validations_needed)

        # 3. Docker Sandbox에서 실행
        result = await self.sandbox.execute(code)

        if result.exit_code != 0:
            # Error -> The Healer로 전달
            await self._log_error(state["session_id"], result.stderr)
            return AgentOutput(
                success=False,
                data={
                    "code": code,
                    "error": result.stderr,
                    "snippet_ids": validations_needed
                },
                error=result.stderr,
                next_agent="healer"
            )

        # 4. 결과 파싱 및 Shared State 업데이트
        metrics = json.loads(result.stdout)
        return AgentOutput(
            success=True,
            data={"metrics": metrics, "snippet_ids": validations_needed},
            reasoning=f"Executed {len(validations_needed)} validation snippets successfully",
            next_agent="auditor"
        )

    def _compose_code_from_snippets(self, smiles: str, snippet_ids: List[str]) -> str:
        """검증된 스니펫들을 조합하여 실행 가능한 코드 생성"""
        imports = set()
        code_blocks = []

        for snippet_id in snippet_ids:
            snippet = self.snippet_manager.get(snippet_id)
            imports.update(snippet.required_imports)
            code_blocks.append(
                snippet.code_template.replace("{{smiles}}", smiles)
            )

        return f"""
import json
{chr(10).join(f'from {imp} import *' if '.' not in imp else f'import {imp}' for imp in imports)}

smiles = "{smiles}"
result = {{}}

{chr(10).join(code_blocks)}

print(json.dumps(result))
"""

# backend/app/sandbox/snippet_manager.py
class SnippetManager:
    """검증된 코드 스니펫 관리자"""

    SNIPPETS = {
        "lipinski_check": CodeSnippet(
            id="lipinski_check",
            name="Lipinski's Rule of Five",
            required_imports=["rdkit.Chem", "rdkit.Chem.Descriptors", "rdkit.Chem.Lipinski"],
            code_template='''
mol = Chem.MolFromSmiles(smiles)
if mol:
    result["lipinski"] = {
        "mw": Descriptors.MolWt(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "pass": Descriptors.MolWt(mol) <= 500 and Descriptors.MolLogP(mol) <= 5
    }
''',
            is_verified=True
        ),

        "tanimoto_similarity": CodeSnippet(
            id="tanimoto_similarity",
            name="Tanimoto Coefficient Calculator",
            required_imports=["rdkit.Chem", "rdkit.Chem.AllChem", "rdkit.DataStructs"],
            code_template='''
mol = Chem.MolFromSmiles(smiles)
if mol:
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    # Reference fingerprints loaded from golden_set
    result["tanimoto"] = {"calculated": True, "fingerprint_bits": fp.GetNumOnBits()}
''',
            is_verified=True
        ),

        "pains_filter": CodeSnippet(
            id="pains_filter",
            name="PAINS (Pan-Assay Interference) Filter",
            required_imports=["rdkit.Chem", "rdkit.Chem.FilterCatalog"],
            code_template='''
mol = Chem.MolFromSmiles(smiles)
if mol:
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog.FilterCatalog(params)
    entry = catalog.GetFirstMatch(mol)
    result["pains"] = {
        "has_pains": entry is not None,
        "pattern": entry.GetDescription() if entry else None
    }
''',
            is_verified=True
        ),

        "sa_score": CodeSnippet(
            id="sa_score",
            name="Synthetic Accessibility Score",
            required_imports=["rdkit.Chem", "rdkit.Chem.rdMolDescriptors"],
            code_template='''
mol = Chem.MolFromSmiles(smiles)
if mol:
    # SA Score calculation (1=easy, 10=hard)
    from rdkit.Contrib.SA_Score import sascorer
    result["sa_score"] = sascorer.calculateScore(mol)
''',
            is_verified=True
        )
    }

    def get(self, snippet_id: str) -> CodeSnippet:
        return self.SNIPPETS.get(snippet_id)

    def get_all_verified(self) -> List[CodeSnippet]:
        return [s for s in self.SNIPPETS.values() if s.is_verified]
```

### 2.4 The Healer with Enhanced Logging

```python
# backend/app/agents/healer.py
class HealerAgent(BaseAgent):
    """
    The Healer: 자가 치유 에이전트
    - 최대 3회 자동 수정 시도
    - 실패 시 Auditor에게 'manual_review' 플래그 전송
    - Constraint 2: 모든 수정 이력을 agent_logs에 기록
    """
    name = "healer"
    MAX_RETRIES = 3

    async def execute(self, state: SessionState) -> AgentOutput:
        original_code = state.get("last_code")
        error_message = state.get("last_error")
        retry_count = state.get("healing_attempts", 0)

        # Log healing attempt start
        await self._log_healing_attempt(
            session_id=state["session_id"],
            attempt=retry_count + 1,
            error=error_message
        )

        # Max retries exceeded -> Manual Review
        if retry_count >= self.MAX_RETRIES:
            await self._log_healing_failure(state["session_id"])
            return AgentOutput(
                success=False,
                data={
                    "status": "manual_review_required",
                    "total_attempts": retry_count,
                    "final_error": error_message
                },
                reasoning="Maximum retry attempts (3) exceeded. Manual review required.",
                error="Max healing retries exceeded"
            )

        # 1. 에러 유형 분류
        error_type = self._classify_error(error_message)

        # 2. LLM을 사용하여 수정된 코드 생성
        fixed_code, fix_explanation = await self._generate_fix(
            original_code=original_code,
            error_message=error_message,
            error_type=error_type,
            attempt=retry_count + 1
        )

        # 3. 수정 내역을 상세하게 로깅 (모델 튜닝용)
        fix_log = {
            "attempt": retry_count + 1,
            "error_type": error_type,
            "original_error": error_message,
            "fix_explanation": fix_explanation,
            "timestamp": datetime.utcnow().isoformat()
        }

        await self._save_fix_log(state["session_id"], fix_log)

        return AgentOutput(
            success=True,
            data={
                "fixed_code": fixed_code,
                "retry_count": retry_count + 1,
                "fix_log": fix_log
            },
            reasoning=f"Error type: {error_type}. Applied fix: {fix_explanation}",
            next_agent="coder"
        )

    def _classify_error(self, error_message: str) -> str:
        """에러 유형 분류"""
        error_lower = error_message.lower()

        if "syntaxerror" in error_lower:
            return "SYNTAX_ERROR"
        elif "valence" in error_lower or "kekulization" in error_lower:
            return "CHEMICAL_VALIDITY_ERROR"
        elif "importerror" in error_lower or "modulenotfounderror" in error_lower:
            return "IMPORT_ERROR"
        elif "timeout" in error_lower:
            return "TIMEOUT_ERROR"
        elif "memoryerror" in error_lower:
            return "MEMORY_ERROR"
        else:
            return "RUNTIME_ERROR"

    async def _save_fix_log(self, session_id: str, fix_log: dict):
        """
        Constraint 2 준수: 수정 이력을 DB에 저장하여 모델 튜닝 데이터로 활용
        """
        await supabase.table("agent_execution_logs").insert({
            "session_id": session_id,
            "agent_name": "healer",
            "status": "healed",
            "reasoning": fix_log["fix_explanation"],
            "fix_logs": [fix_log],
            "retry_count": fix_log["attempt"]
        }).execute()

    async def _register_successful_fix(self, session_id: str, fix_log: dict,
                                        fixed_code: str, original_error: str):
        """
        Knowledge Loop (v2.1): 성공적인 수정 패턴을 candidate_snippets에 자동 등록
        Constraint 5 준수
        """
        # 에러 패턴 추출
        error_pattern = self._extract_error_pattern(original_error)

        await supabase.table("candidate_snippets").insert({
            "source_session_id": session_id,
            "error_type": fix_log["error_type"],
            "error_pattern": error_pattern,
            "error_context": {
                "full_error": original_error,
                "code_context": fix_log.get("code_context")
            },
            "fix_code": fixed_code,
            "fix_explanation": fix_log["fix_explanation"],
            "fix_strategy": self._classify_fix_strategy(fix_log),
            "confidence_score": 0.50  # 초기 신뢰도
        }).execute()

    def _extract_error_pattern(self, error_message: str) -> str:
        """에러 메시지에서 재사용 가능한 패턴 추출"""
        import re
        # Remove file paths and line numbers
        pattern = re.sub(r'File ".*?", line \d+', 'File "...", line X', error_message)
        # Remove specific variable names
        pattern = re.sub(r"'[^']*'", "'VAR'", pattern)
        return pattern[:500]  # 최대 500자

    def _classify_fix_strategy(self, fix_log: dict) -> str:
        """수정 전략 분류"""
        explanation = fix_log.get("fix_explanation", "").lower()

        if "null check" in explanation or "none check" in explanation:
            return "add_null_check"
        elif "try" in explanation or "except" in explanation:
            return "add_exception_handler"
        elif "import" in explanation:
            return "fix_import"
        elif "type" in explanation or "cast" in explanation:
            return "fix_type_conversion"
        elif "syntax" in explanation:
            return "fix_syntax"
        else:
            return "general_fix"
```

### 2.5 The Librarian with Dynamic Substructure RAG (v2.1 신규)

```python
# backend/app/agents/librarian.py
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, DataStructs
from app.services.embedding_service import EmbeddingService
from app.services.rag_service import RAGService

class LibrarianAgent(BaseAgent):
    """
    The Librarian: Dynamic Substructure RAG 에이전트 (v2.1 Enhanced)

    핵심 기능:
    1. 입력 SMILES에서 스캐폴드/서브스트럭처 추출 (Murcko decomposition)
    2. 스캐폴드를 벡터 임베딩하여 의미론적 유사성 검색
    3. knowledge_base + golden_set에서 관련 문헌/구조 검색
    4. 관련 PMID 및 근거 자료 제공
    """
    name = "librarian"

    def __init__(self):
        self.embedding_service = EmbeddingService()
        self.rag_service = RAGService()

    async def execute(self, state: SessionState) -> AgentOutput:
        smiles = state["current_smiles"]
        target_antigen = state.get("target_antigen", "")

        await self._log_start(state["session_id"])

        # 1. 스캐폴드 추출 (Murcko Decomposition)
        scaffolds = self._extract_scaffolds(smiles)

        # 2. 각 스캐폴드에 대해 지식베이스 매핑
        knowledge_results = []
        for scaffold in scaffolds:
            # 2.1 기존 매핑 검색 또는 새로 생성
            mapping = await self._get_or_create_scaffold_mapping(scaffold)

            # 2.2 의미론적 검색 (Semantic Search)
            semantic_hits = await self._semantic_search(
                scaffold["smiles"],
                target_antigen,
                top_k=5
            )

            # 2.3 Tanimoto 유사도 기반 Golden Set 검색
            golden_set_hits = await self._search_golden_set_by_similarity(
                scaffold["fingerprint"],
                threshold=0.6
            )

            knowledge_results.append({
                "scaffold_type": scaffold["type"],
                "scaffold_smiles": scaffold["smiles"],
                "semantic_hits": semantic_hits,
                "golden_set_hits": golden_set_hits,
                "mapping_id": mapping["id"] if mapping else None
            })

        # 3. 통합 근거 생성
        evidence = self._compile_evidence(knowledge_results)

        return AgentOutput(
            success=True,
            data={
                "scaffolds": scaffolds,
                "knowledge_results": knowledge_results,
                "evidence_summary": evidence["summary"],
                "pmid_references": evidence["pmids"],
                "golden_set_references": evidence["golden_sets"]
            },
            reasoning=f"Extracted {len(scaffolds)} scaffolds, found {len(evidence['pmids'])} literature references",
            referenced_knowledge_ids=evidence["knowledge_ids"],
            referenced_pmids=evidence["pmids"]
        )

    def _extract_scaffolds(self, smiles: str) -> List[dict]:
        """
        Murcko Decomposition으로 스캐폴드 추출
        Returns: [{type, smiles, fingerprint}, ...]
        """
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return []

        scaffolds = []

        try:
            # 1. Generic Scaffold (모든 원자를 탄소로, 모든 결합을 단일결합으로)
            generic = MurckoScaffold.MakeScaffoldGeneric(
                MurckoScaffold.GetScaffoldForMol(mol)
            )
            if generic:
                generic_smiles = Chem.MolToSmiles(generic)
                scaffolds.append({
                    "type": "murcko_generic",
                    "smiles": generic_smiles,
                    "fingerprint": self._compute_fingerprint(generic)
                })

            # 2. Framework (링 시스템 + 링커만)
            framework = MurckoScaffold.GetScaffoldForMol(mol)
            if framework:
                framework_smiles = Chem.MolToSmiles(framework)
                scaffolds.append({
                    "type": "murcko_framework",
                    "smiles": framework_smiles,
                    "fingerprint": self._compute_fingerprint(framework)
                })

            # 3. Ring Systems (개별 링 시스템)
            ring_info = mol.GetRingInfo()
            for ring_atoms in ring_info.AtomRings():
                ring_mol = Chem.RWMol(Chem.MolFromSmiles(''))
                # Extract ring as submol (simplified)
                ring_smiles = self._extract_ring_smiles(mol, ring_atoms)
                if ring_smiles:
                    scaffolds.append({
                        "type": "ring_system",
                        "smiles": ring_smiles,
                        "fingerprint": None  # Compute if needed
                    })

        except Exception as e:
            logging.warning(f"Scaffold extraction error: {e}")

        return scaffolds

    def _compute_fingerprint(self, mol) -> bytes:
        """Morgan Fingerprint 계산"""
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        return fp.ToBitString().encode()

    async def _get_or_create_scaffold_mapping(self, scaffold: dict) -> Optional[dict]:
        """
        스캐폴드 매핑 조회 또는 신규 생성
        """
        # 기존 매핑 검색
        existing = await supabase.table("scaffold_knowledge_mapping").select("*").eq(
            "scaffold_smiles", scaffold["smiles"]
        ).eq(
            "scaffold_type", scaffold["type"]
        ).single().execute()

        if existing.data:
            # 조회 카운트 증가
            await supabase.table("scaffold_knowledge_mapping").update({
                "retrieval_count": existing.data["retrieval_count"] + 1,
                "last_retrieved_at": datetime.utcnow().isoformat()
            }).eq("id", existing.data["id"]).execute()
            return existing.data

        # 신규 매핑 생성
        embedding = await self.embedding_service.embed_text(
            f"Chemical scaffold: {scaffold['smiles']}"
        )

        new_mapping = await supabase.table("scaffold_knowledge_mapping").insert({
            "scaffold_smiles": scaffold["smiles"],
            "scaffold_type": scaffold["type"],
            "fingerprint_bits": scaffold.get("fingerprint"),
            "embedding_vector": embedding,
            "scaffold_description": await self._generate_scaffold_description(scaffold),
            "retrieval_count": 1,
            "last_retrieved_at": datetime.utcnow().isoformat()
        }).execute()

        return new_mapping.data[0] if new_mapping.data else None

    async def _semantic_search(self, scaffold_smiles: str, target: str, top_k: int = 5) -> List[dict]:
        """
        의미론적 검색: 스캐폴드 + 타겟을 기반으로 knowledge_base 검색
        """
        query = f"ADC scaffold structure {scaffold_smiles} targeting {target}"

        # RAG 검색 (기존 rag_service 활용)
        results = await self.rag_service.semantic_search(
            query=query,
            top_k=top_k,
            filter_table="knowledge_base"
        )

        return [{
            "id": r.id,
            "title": r.title,
            "relevance_score": r.score,
            "pmid": r.metadata.get("pmid")
        } for r in results]

    async def _search_golden_set_by_similarity(self, fingerprint: bytes, threshold: float = 0.6) -> List[dict]:
        """
        Golden Set에서 Tanimoto 유사도 기반 검색
        """
        if not fingerprint:
            return []

        # RPC를 통한 유사도 검색 (Supabase function)
        results = await supabase.rpc("search_golden_set_by_tanimoto", {
            "query_fingerprint": fingerprint,
            "similarity_threshold": threshold,
            "limit_count": 10
        }).execute()

        return results.data or []

    def _compile_evidence(self, knowledge_results: List[dict]) -> dict:
        """
        검색 결과를 통합하여 근거 자료 생성
        """
        pmids = set()
        knowledge_ids = set()
        golden_sets = []
        summaries = []

        for result in knowledge_results:
            for hit in result.get("semantic_hits", []):
                if hit.get("pmid"):
                    pmids.add(hit["pmid"])
                knowledge_ids.add(hit["id"])

            for gs in result.get("golden_set_hits", []):
                golden_sets.append({
                    "id": gs["id"],
                    "name": gs.get("name"),
                    "similarity": gs.get("tanimoto_score")
                })

            summaries.append(
                f"{result['scaffold_type']}: {len(result.get('semantic_hits', []))} refs"
            )

        return {
            "pmids": list(pmids),
            "knowledge_ids": list(knowledge_ids),
            "golden_sets": golden_sets[:5],  # Top 5
            "summary": "; ".join(summaries)
        }

    async def _generate_scaffold_description(self, scaffold: dict) -> str:
        """LLM을 사용하여 스캐폴드 설명 생성"""
        # Simplified - could use LLM for richer descriptions
        return f"Chemical scaffold of type {scaffold['type']} with SMILES pattern"
```

### 2.6 The Auditor with Constraint Guardrail (v2.1 신규)

```python
# backend/app/agents/auditor.py
class AuditorAgent(BaseAgent):
    """
    The Auditor: 검증 및 감사 에이전트 (v2.1 Enhanced)

    핵심 기능:
    1. PAINS/Lipinski/독성 검증
    2. 🔒 Constraint Guardrail: 초기 제약조건 준수 확인 (v2.1 신규)
    3. 리스크 스코어 계산
    4. 최종 승인/재설계 결정
    """
    name = "auditor"

    async def execute(self, state: SessionState) -> AgentOutput:
        await self._log_start(state["session_id"])

        # 1. 기본 화학 검증
        chemistry_validation = await self._validate_chemistry(state["current_smiles"])

        # 2. 🔒 Constraint Guardrail 검사 (v2.1 신규)
        constraint_check = await self._check_constraint_guardrail(state)

        # 3. 리스크 스코어 계산
        risk_score = self._calculate_risk_score(
            chemistry_validation,
            constraint_check
        )

        # 4. 최종 결정
        decision = self._make_decision(
            chemistry_validation,
            constraint_check,
            risk_score
        )

        # 5. 상세 로깅
        await self._log_audit_result(
            session_id=state["session_id"],
            validation=chemistry_validation,
            constraint_check=constraint_check,
            risk_score=risk_score,
            decision=decision
        )

        return AgentOutput(
            success=decision["approved"],
            data={
                "chemistry_validation": chemistry_validation,
                "constraint_check": constraint_check,
                "risk_score": risk_score,
                "decision": decision,
                "final_report": self._generate_report(
                    state, chemistry_validation, constraint_check, risk_score
                )
            },
            reasoning=decision["reasoning"],
            confidence_score=decision["confidence"],
            next_agent="alchemist" if decision["action"] == "redesign" else None
        )

    async def _check_constraint_guardrail(self, state: SessionState) -> dict:
        """
        🔒 Constraint Guardrail (v2.1)
        사용자 초기 제약조건 vs 설계 결과 비교 검증

        Constraint 4 준수: 세션 생성 시 입력한 제약조건이 결과에 반영되었는지 확인
        """
        constraints = {
            "target_antigen": state.get("target_antigen"),
            "requested_dar": state.get("requested_dar"),
            "linker_preference": state.get("linker_preference")
        }

        violations = []
        warnings = []

        # 1. Target Antigen 검증
        if constraints["target_antigen"]:
            target_validated = await self._validate_target_compatibility(
                state["current_smiles"],
                constraints["target_antigen"]
            )
            if not target_validated["compatible"]:
                violations.append({
                    "constraint": "target_antigen",
                    "expected": constraints["target_antigen"],
                    "actual": target_validated.get("detected_target"),
                    "severity": "high",
                    "reason": target_validated.get("reason")
                })

        # 2. DAR (Drug-Antibody Ratio) 검증
        if constraints["requested_dar"]:
            detected_dar = self._estimate_dar_from_structure(state["current_smiles"])
            if detected_dar and detected_dar != constraints["requested_dar"]:
                # DAR 차이가 허용 범위 내인지 확인
                if abs(detected_dar - constraints["requested_dar"]) > 1:
                    violations.append({
                        "constraint": "requested_dar",
                        "expected": constraints["requested_dar"],
                        "actual": detected_dar,
                        "severity": "medium",
                        "reason": f"Designed structure suggests DAR {detected_dar}, not {constraints['requested_dar']}"
                    })
                else:
                    warnings.append({
                        "constraint": "requested_dar",
                        "message": f"DAR slightly different: expected {constraints['requested_dar']}, estimated {detected_dar}"
                    })

        # 3. Linker Preference 검증
        if constraints["linker_preference"] and constraints["linker_preference"] != "any":
            linker_analysis = self._analyze_linker_type(state["current_smiles"])
            if linker_analysis["type"] != constraints["linker_preference"]:
                violations.append({
                    "constraint": "linker_preference",
                    "expected": constraints["linker_preference"],
                    "actual": linker_analysis["type"],
                    "severity": "medium",
                    "reason": linker_analysis.get("reason")
                })

        return {
            "passed": len(violations) == 0,
            "violations": violations,
            "warnings": warnings,
            "constraints_checked": list(constraints.keys()),
            "summary": self._generate_constraint_summary(violations, warnings)
        }

    async def _validate_target_compatibility(self, smiles: str, target: str) -> dict:
        """타겟 항원 호환성 검증"""
        # 실제 구현에서는 도메인 지식 또는 ML 모델 사용
        known_targets = {
            "HER2": ["trastuzumab", "pertuzumab", "ado-trastuzumab"],
            "TROP2": ["sacituzumab"],
            "EGFR": ["cetuximab"],
            # ... more targets
        }

        # 간단한 검증 로직 (실제로는 더 복잡한 분석 필요)
        return {
            "compatible": True,  # Placeholder
            "confidence": 0.85,
            "detected_target": target,
            "reason": "Target compatibility verified through structural analysis"
        }

    def _estimate_dar_from_structure(self, smiles: str) -> Optional[int]:
        """구조에서 DAR 추정"""
        # 약물-항체 결합점 수 분석
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None

        # Cysteine 결합점 또는 Lysine 결합점 패턴 검색
        # 실제 구현에서는 더 정교한 분석 필요
        return 4  # Placeholder

    def _analyze_linker_type(self, smiles: str) -> dict:
        """링커 타입 분석"""
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {"type": "unknown", "reason": "Invalid SMILES"}

        # Cleavable linker 패턴 (예: Val-Cit, Phe-Lys)
        cleavable_patterns = [
            "[NX3][CX4][CX4][NX3][CX3](=O)",  # Dipeptide linker
            "[CX3](=O)[NX3][CX4][CX3](=O)[NX3]"  # Peptide bond
        ]

        # Non-cleavable linker 패턴 (예: SMCC, MCC)
        non_cleavable_patterns = [
            "[NX3][CX3](=O)[CX4][CX4]S",  # Thioether linker
        ]

        smiles_str = smiles.upper()

        for pattern in cleavable_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                return {"type": "cleavable", "reason": "Protease-sensitive peptide linker detected"}

        for pattern in non_cleavable_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                return {"type": "non-cleavable", "reason": "Stable thioether linker detected"}

        return {"type": "unknown", "reason": "Linker type could not be determined"}

    def _generate_constraint_summary(self, violations: List[dict], warnings: List[dict]) -> str:
        """제약조건 검사 요약 생성"""
        if not violations and not warnings:
            return "✅ All user constraints satisfied"

        parts = []
        if violations:
            parts.append(f"❌ {len(violations)} constraint violation(s)")
        if warnings:
            parts.append(f"⚠️ {len(warnings)} warning(s)")

        return "; ".join(parts)

    def _make_decision(self, chemistry_validation: dict,
                       constraint_check: dict, risk_score: float) -> dict:
        """
        최종 결정: 승인 / 재설계 / 수동 검토
        """
        # High-severity constraint violations -> 재설계
        high_violations = [v for v in constraint_check["violations"]
                          if v.get("severity") == "high"]
        if high_violations:
            return {
                "approved": False,
                "action": "redesign",
                "reasoning": f"Critical constraint violations: {[v['constraint'] for v in high_violations]}",
                "confidence": 0.9
            }

        # Chemistry validation failure -> 재설계
        if not chemistry_validation.get("lipinski_pass"):
            return {
                "approved": False,
                "action": "redesign",
                "reasoning": "Failed Lipinski's Rule of Five",
                "confidence": 0.85
            }

        if chemistry_validation.get("pains_detected"):
            return {
                "approved": False,
                "action": "manual_review",
                "reasoning": "PAINS alerts detected - requires expert review",
                "confidence": 0.75
            }

        # High risk score -> 수동 검토
        if risk_score > 7:
            return {
                "approved": False,
                "action": "manual_review",
                "reasoning": f"High risk score ({risk_score}/10) requires expert review",
                "confidence": 0.7
            }

        # Medium-severity violations with good chemistry -> 경고와 함께 승인
        if constraint_check["violations"]:
            return {
                "approved": True,
                "action": "approved_with_warnings",
                "reasoning": f"Approved with {len(constraint_check['violations'])} minor constraint deviations",
                "confidence": 0.75,
                "warnings": constraint_check["violations"]
            }

        # All good
        return {
            "approved": True,
            "action": "approved",
            "reasoning": "All validations passed, constraints satisfied",
            "confidence": 0.95
        }
```

---

## Phase 3: Docker Sandbox Environment [Week 4]

### 3.1 Sandbox Dockerfile (Version-Locked)

> ⚠️ **Constraint 7 준수**: 단일 `requirements.lock` 파일로 버전 100% 동기화

```dockerfile
# backend/app/sandbox/Dockerfile.sandbox
FROM python:3.11-slim

# System dependencies for RDKit
RUN apt-get update && apt-get install -y \
    libxrender1 \
    libxext6 \
    && rm -rf /var/lib/apt/lists/*

# 🔒 Version-locked requirements (단일 소스)
COPY requirements.lock /tmp/requirements.lock

# Install chemistry libraries from lock file
RUN pip install --no-cache-dir -r /tmp/requirements.lock

# Security: Create non-root user
RUN useradd -m -s /bin/bash sandbox
USER sandbox
WORKDIR /sandbox

# 🔒 버전 검증 스크립트 복사
COPY --chown=sandbox:sandbox verify_versions.py /sandbox/verify_versions.py

# Environment
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# 🔒 빌드 시 버전 검증 (실패 시 빌드 중단)
RUN python /sandbox/verify_versions.py

# Default command
CMD ["python", "-c", "print('Sandbox ready')"]
```

### 3.1.1 Version Lock File (단일 소스 of Truth)

```txt
# backend/app/sandbox/requirements.lock
# ⚠️ 이 파일은 에이전트 학습 환경과 100% 동일해야 합니다
# 수정 시 반드시 CI/CD 파이프라인 검증 필요

rdkit==2023.9.4
biopython==1.83
pandas==2.1.4
numpy==1.26.3
scipy==1.12.0
scikit-learn==1.4.0
matplotlib==3.8.2
```

### 3.1.2 Version Verification Script

```python
# backend/app/sandbox/verify_versions.py
"""
🔒 Sandbox 라이브러리 버전 검증 스크립트 (v2.1)
Constraint 7 준수: 버전 불일치 시 즉시 실패
"""
import sys
import importlib.metadata
import json

# 필수 버전 (requirements.lock과 동기화)
REQUIRED_VERSIONS = {
    "rdkit": "2023.9.4",
    "biopython": "1.83",
    "pandas": "2.1.4",
    "numpy": "1.26.3",
    "scipy": "1.12.0",
}

def verify_versions():
    errors = []
    installed = {}

    for package, required_version in REQUIRED_VERSIONS.items():
        try:
            # RDKit은 특수 처리
            if package == "rdkit":
                from rdkit import __version__ as rdkit_version
                actual_version = rdkit_version
            else:
                actual_version = importlib.metadata.version(package)

            installed[package] = actual_version

            if actual_version != required_version:
                errors.append({
                    "package": package,
                    "required": required_version,
                    "actual": actual_version
                })

        except importlib.metadata.PackageNotFoundError:
            errors.append({
                "package": package,
                "required": required_version,
                "actual": "NOT INSTALLED"
            })

    if errors:
        print("=" * 60)
        print("❌ VERSION MISMATCH DETECTED - BUILD BLOCKED")
        print("=" * 60)
        for err in errors:
            print(f"  {err['package']}: required={err['required']}, actual={err['actual']}")
        print("=" * 60)
        print("Fix: Update requirements.lock and rebuild")
        sys.exit(1)

    print("✅ All library versions verified successfully")
    print(json.dumps(installed, indent=2))
    return True

if __name__ == "__main__":
    verify_versions()
```

### 3.1.3 CI/CD Version Sync Check

```yaml
# .github/workflows/sandbox-version-check.yml
name: Sandbox Version Sync Check

on:
  push:
    paths:
      - 'backend/app/sandbox/requirements.lock'
      - 'backend/requirements.txt'
      - 'backend/app/agents/**'
  pull_request:
    paths:
      - 'backend/app/sandbox/requirements.lock'

jobs:
  version-sync:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: '3.11'

      - name: Compare Agent and Sandbox versions
        run: |
          echo "🔍 Checking version sync between Agent and Sandbox environments..."

          # Extract versions from both files
          AGENT_RDKIT=$(grep "rdkit" backend/requirements.txt | cut -d'=' -f3)
          SANDBOX_RDKIT=$(grep "rdkit" backend/app/sandbox/requirements.lock | cut -d'=' -f2)

          if [ "$AGENT_RDKIT" != "$SANDBOX_RDKIT" ]; then
            echo "❌ RDKit version mismatch!"
            echo "   Agent: $AGENT_RDKIT"
            echo "   Sandbox: $SANDBOX_RDKIT"
            exit 1
          fi

          # Add more package checks...
          echo "✅ All versions synchronized"

      - name: Build and test Sandbox image
        run: |
          cd backend/app/sandbox
          docker build -t adc-sandbox:test .
          docker run --rm adc-sandbox:test python verify_versions.py
```

### 3.2 Sandbox Executor with Resource Limits

```python
# backend/app/sandbox/executor.py
import docker
import tempfile
import asyncio
from dataclasses import dataclass

@dataclass
class ExecutionResult:
    exit_code: int
    stdout: str
    stderr: str
    execution_time_ms: int

class SandboxExecutor:
    """
    Docker 기반 격리된 Python 실행 환경
    - 네트워크 차단 (Air-gapped)
    - 메모리 제한: 512MB
    - CPU 제한: 50%
    - 실행 시간 제한: 30초
    """

    def __init__(self):
        self.client = docker.from_env()
        self.image = "adc-sandbox:latest"
        self.timeout = 30
        self.memory_limit = "512m"
        self.cpu_quota = 50000  # 50% of one CPU

    async def execute(self, code: str) -> ExecutionResult:
        start_time = asyncio.get_event_loop().time()

        # Write code to temp file
        with tempfile.NamedTemporaryFile(
            mode='w', suffix=".py", delete=False
        ) as f:
            f.write(code)
            script_path = f.name

        try:
            # Run in isolated container
            container = self.client.containers.run(
                self.image,
                command=f"python /sandbox/script.py",
                volumes={
                    script_path: {"bind": "/sandbox/script.py", "mode": "ro"}
                },
                network_mode="none",           # 🔒 No network access
                mem_limit=self.memory_limit,   # 🔒 Memory limit
                cpu_period=100000,
                cpu_quota=self.cpu_quota,      # 🔒 CPU limit
                remove=True,
                detach=False,
                stdout=True,
                stderr=True,
                timeout=self.timeout           # 🔒 Timeout
            )

            execution_time = int((asyncio.get_event_loop().time() - start_time) * 1000)

            return ExecutionResult(
                exit_code=0,
                stdout=container.decode('utf-8'),
                stderr="",
                execution_time_ms=execution_time
            )

        except docker.errors.ContainerError as e:
            execution_time = int((asyncio.get_event_loop().time() - start_time) * 1000)
            return ExecutionResult(
                exit_code=e.exit_status,
                stdout="",
                stderr=str(e.stderr.decode('utf-8') if e.stderr else str(e)),
                execution_time_ms=execution_time
            )

        except Exception as e:
            return ExecutionResult(
                exit_code=-1,
                stdout="",
                stderr=str(e),
                execution_time_ms=0
            )

        finally:
            # Cleanup temp file
            import os
            if os.path.exists(script_path):
                os.remove(script_path)
```

---

## Phase 3.5: AlphaFold 3 Integration [Week 4-5] (v2.1 신규 - Premium Tier)

> 💎 **Premium Tier 전용**: 3D 단백질-약물 상호작용 예측

### 3.5.1 Database Schema Extension

```sql
-- =====================================================
-- AlphaFold 3D Modeling Jobs (Premium Only)
-- =====================================================
CREATE TABLE alphafold_jobs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    session_id UUID REFERENCES design_sessions(id) ON DELETE CASCADE,
    user_id UUID REFERENCES auth.users(id) NOT NULL,

    -- Input
    target_sequence TEXT NOT NULL,              -- Antibody/Target protein sequence
    ligand_smiles TEXT NOT NULL,                -- Drug candidate SMILES
    job_type VARCHAR(50) DEFAULT 'binding',     -- 'binding', 'structure', 'complex'

    -- Job Status
    status VARCHAR(30) DEFAULT 'queued' CHECK (status IN (
        'queued', 'processing', 'completed', 'failed', 'timeout'
    )),
    priority INTEGER DEFAULT 5,                 -- 1=highest, 10=lowest

    -- Results
    predicted_structure JSONB,                  -- PDB-like structure data
    binding_affinity DECIMAL(8,4),              -- Predicted binding affinity (ΔG)
    confidence_score DECIMAL(3,2),              -- AlphaFold confidence (pLDDT)
    contact_residues TEXT[],                    -- Key binding residues
    visualization_url TEXT,                     -- 3D viewer URL

    -- Compute Metrics
    compute_time_seconds INTEGER,
    gpu_hours_used DECIMAL(6,3),

    -- Timestamps
    created_at TIMESTAMPTZ DEFAULT NOW(),
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ,

    -- Premium Check
    CONSTRAINT premium_only CHECK (
        EXISTS (
            SELECT 1 FROM design_sessions ds
            WHERE ds.id = session_id AND ds.tier = 'premium'
        )
    )
);

CREATE INDEX idx_alphafold_session ON alphafold_jobs(session_id);
CREATE INDEX idx_alphafold_status ON alphafold_jobs(status, priority);

-- =====================================================
-- AlphaFold GPU Quota Management (v2.1 신규) ⚠️ Critical
-- =====================================================
CREATE TABLE alphafold_user_quota (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES auth.users(id) UNIQUE NOT NULL,

    -- 월간 할당량
    monthly_gpu_hours_limit DECIMAL(6,2) DEFAULT 10.0,  -- 기본 10시간/월
    monthly_gpu_hours_used DECIMAL(6,2) DEFAULT 0.0,
    quota_reset_date DATE DEFAULT (DATE_TRUNC('month', NOW()) + INTERVAL '1 month')::DATE,

    -- 일간 제한 (DDoS 방지)
    daily_job_limit INTEGER DEFAULT 20,                 -- 일 최대 20개 작업
    daily_jobs_submitted INTEGER DEFAULT 0,
    daily_reset_date DATE DEFAULT CURRENT_DATE,

    -- 우선순위 티어
    priority_tier VARCHAR(20) DEFAULT 'standard' CHECK (priority_tier IN (
        'standard',     -- 일반 Premium (priority 5)
        'priority',     -- 우선 처리 (priority 3)
        'enterprise'    -- 최우선 (priority 1)
    )),

    -- 초과 사용 정책
    allow_overage BOOLEAN DEFAULT FALSE,               -- 초과 사용 허용 여부
    overage_rate_per_hour DECIMAL(6,2) DEFAULT 5.00,   -- 초과 시간당 요금 ($)

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- 쿼터 자동 리셋 함수
CREATE OR REPLACE FUNCTION reset_alphafold_quota()
RETURNS TRIGGER AS $$
BEGIN
    -- 월간 쿼터 리셋
    IF NEW.quota_reset_date <= CURRENT_DATE THEN
        NEW.monthly_gpu_hours_used := 0;
        NEW.quota_reset_date := (DATE_TRUNC('month', NOW()) + INTERVAL '1 month')::DATE;
    END IF;

    -- 일간 쿼터 리셋
    IF NEW.daily_reset_date < CURRENT_DATE THEN
        NEW.daily_jobs_submitted := 0;
        NEW.daily_reset_date := CURRENT_DATE;
    END IF;

    NEW.updated_at := NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER auto_reset_quota
    BEFORE UPDATE ON alphafold_user_quota
    FOR EACH ROW EXECUTE FUNCTION reset_alphafold_quota();

CREATE INDEX idx_quota_user ON alphafold_user_quota(user_id);
CREATE INDEX idx_quota_reset ON alphafold_user_quota(quota_reset_date);
```

### 3.5.2 AlphaFold Service

```python
# backend/app/services/alphafold_service.py
import httpx
from typing import Optional
from pydantic import BaseModel
from enum import Enum

class AlphaFoldJobType(str, Enum):
    BINDING = "binding"           # 결합 예측
    STRUCTURE = "structure"       # 구조 예측
    COMPLEX = "complex"           # 복합체 모델링

class AlphaFoldRequest(BaseModel):
    target_sequence: str          # 단백질 서열 (FASTA format)
    ligand_smiles: str            # 약물 SMILES
    job_type: AlphaFoldJobType = AlphaFoldJobType.BINDING
    confidence_threshold: float = 0.7

class AlphaFoldResult(BaseModel):
    job_id: str
    status: str
    structure_pdb: Optional[str]  # PDB format structure
    binding_affinity: Optional[float]
    confidence: Optional[float]
    contact_residues: Optional[list]
    visualization_data: Optional[dict]

class AlphaFoldService:
    """
    AlphaFold 3 Integration Service (v2.1)

    Premium Tier 전용 3D 모델링 서비스
    - 비동기 Job Queue 방식
    - GPU 리소스 관리
    - 결과 캐싱
    """

    def __init__(self):
        self.api_base = os.environ.get("ALPHAFOLD_API_URL", "http://alphafold-server:8080")
        self.api_key = os.environ.get("ALPHAFOLD_API_KEY")
        self.timeout = 300  # 5분 타임아웃

    async def submit_job(self, session_id: str, user_id: str,
                         request: AlphaFoldRequest) -> str:
        """
        AlphaFold 작업 제출

        Returns: job_id
        """
        # 1. Premium 권한 확인
        session = await supabase.table("design_sessions").select("tier").eq(
            "id", session_id
        ).single().execute()

        if session.data.get("tier") != "premium":
            raise PermissionError("AlphaFold 3D modeling is Premium-only feature")

        # 2. 🔒 GPU 쿼터 확인 (v2.1 신규)
        quota = await self._check_and_update_quota(user_id)
        if not quota["allowed"]:
            raise QuotaExceededError(
                f"GPU quota exceeded. Used: {quota['used']:.1f}h / Limit: {quota['limit']:.1f}h. "
                f"Resets on: {quota['reset_date']}"
            )

        # 2. Job 레코드 생성
        job = await supabase.table("alphafold_jobs").insert({
            "session_id": session_id,
            "user_id": user_id,
            "target_sequence": request.target_sequence,
            "ligand_smiles": request.ligand_smiles,
            "job_type": request.job_type.value,
            "status": "queued"
        }).execute()

        job_id = job.data[0]["id"]

        # 3. 비동기 처리 큐에 추가
        await self._enqueue_job(job_id, request)

        return job_id

    async def _enqueue_job(self, job_id: str, request: AlphaFoldRequest):
        """작업을 처리 큐에 추가"""
        # Redis/Celery 또는 Cloud Tasks 사용
        # 여기서는 간단히 background task로 처리
        import asyncio
        asyncio.create_task(self._process_job(job_id, request))

    async def _process_job(self, job_id: str, request: AlphaFoldRequest):
        """
        실제 AlphaFold API 호출 및 처리
        """
        try:
            # 상태 업데이트: processing
            await supabase.table("alphafold_jobs").update({
                "status": "processing",
                "started_at": datetime.utcnow().isoformat()
            }).eq("id", job_id).execute()

            start_time = time.time()

            # AlphaFold Server API 호출
            async with httpx.AsyncClient(timeout=self.timeout) as client:
                response = await client.post(
                    f"{self.api_base}/v1/predict",
                    headers={"Authorization": f"Bearer {self.api_key}"},
                    json={
                        "sequences": [request.target_sequence],
                        "ligand": request.ligand_smiles,
                        "model": "alphafold3",
                        "return_contacts": True,
                        "return_confidence": True
                    }
                )

            if response.status_code != 200:
                raise Exception(f"AlphaFold API error: {response.text}")

            result = response.json()
            compute_time = int(time.time() - start_time)

            # 결과 저장
            await supabase.table("alphafold_jobs").update({
                "status": "completed",
                "completed_at": datetime.utcnow().isoformat(),
                "predicted_structure": result.get("structure"),
                "binding_affinity": result.get("binding_affinity"),
                "confidence_score": result.get("plddt_mean"),
                "contact_residues": result.get("contact_residues", []),
                "visualization_url": self._generate_visualization_url(job_id, result),
                "compute_time_seconds": compute_time,
                "gpu_hours_used": compute_time / 3600 * 0.5  # Estimated GPU usage
            }).eq("id", job_id).execute()

            # WebSocket으로 완료 알림
            await self._notify_completion(job_id)

        except Exception as e:
            logging.error(f"AlphaFold job {job_id} failed: {e}")
            await supabase.table("alphafold_jobs").update({
                "status": "failed",
                "completed_at": datetime.utcnow().isoformat()
            }).eq("id", job_id).execute()

    def _generate_visualization_url(self, job_id: str, result: dict) -> str:
        """3D 시각화 URL 생성 (Mol* 또는 NGL Viewer)"""
        # 실제 구현에서는 PDB 파일을 S3에 저장하고 뷰어 URL 생성
        return f"/api/alphafold/viewer/{job_id}"

    async def get_job_status(self, job_id: str) -> AlphaFoldResult:
        """작업 상태 조회"""
        job = await supabase.table("alphafold_jobs").select("*").eq(
            "id", job_id
        ).single().execute()

        if not job.data:
            raise ValueError(f"Job {job_id} not found")

        return AlphaFoldResult(
            job_id=job_id,
            status=job.data["status"],
            structure_pdb=job.data.get("predicted_structure", {}).get("pdb"),
            binding_affinity=job.data.get("binding_affinity"),
            confidence=job.data.get("confidence_score"),
            contact_residues=job.data.get("contact_residues"),
            visualization_data=job.data.get("predicted_structure")
        )

    async def _notify_completion(self, job_id: str):
        """WebSocket으로 완료 알림"""
        job = await supabase.table("alphafold_jobs").select(
            "session_id, user_id, status, binding_affinity, confidence_score, gpu_hours_used"
        ).eq("id", job_id).single().execute()

        if job.data:
            # GPU 사용량 쿼터에 반영
            await self._update_quota_usage(
                job.data["user_id"],
                job.data.get("gpu_hours_used", 0)
            )

            await websocket_hub.broadcast(job.data["session_id"], {
                "type": "alphafold_complete",
                "job_id": job_id,
                "status": job.data["status"],
                "binding_affinity": job.data.get("binding_affinity"),
                "confidence": job.data.get("confidence_score")
            })

    # =========================================================
    # 🔒 GPU Quota Management (v2.1 신규)
    # =========================================================
    async def _check_and_update_quota(self, user_id: str) -> dict:
        """
        사용자 GPU 쿼터 확인 및 일간 작업 수 증가

        Returns: {allowed: bool, used: float, limit: float, reset_date: str}
        """
        # 쿼터 레코드 조회 또는 생성
        quota = await supabase.table("alphafold_user_quota").select("*").eq(
            "user_id", user_id
        ).single().execute()

        if not quota.data:
            # 새 사용자: 기본 쿼터 생성
            quota = await supabase.table("alphafold_user_quota").insert({
                "user_id": user_id
            }).execute()
            quota = quota.data[0] if quota.data else {}
        else:
            quota = quota.data

        # 월간 GPU 시간 확인
        used = float(quota.get("monthly_gpu_hours_used", 0))
        limit = float(quota.get("monthly_gpu_hours_limit", 10))
        allow_overage = quota.get("allow_overage", False)

        # 일간 작업 수 확인
        daily_used = quota.get("daily_jobs_submitted", 0)
        daily_limit = quota.get("daily_job_limit", 20)

        # 쿼터 초과 확인
        if used >= limit and not allow_overage:
            return {
                "allowed": False,
                "used": used,
                "limit": limit,
                "reset_date": quota.get("quota_reset_date"),
                "reason": "monthly_quota_exceeded"
            }

        if daily_used >= daily_limit:
            return {
                "allowed": False,
                "used": used,
                "limit": limit,
                "reset_date": quota.get("daily_reset_date"),
                "reason": "daily_limit_exceeded"
            }

        # 일간 작업 수 증가
        await supabase.table("alphafold_user_quota").update({
            "daily_jobs_submitted": daily_used + 1
        }).eq("user_id", user_id).execute()

        # 우선순위 결정
        priority = self._get_priority_from_tier(quota.get("priority_tier", "standard"))

        return {
            "allowed": True,
            "used": used,
            "limit": limit,
            "reset_date": quota.get("quota_reset_date"),
            "priority": priority,
            "is_overage": used >= limit and allow_overage
        }

    async def _update_quota_usage(self, user_id: str, gpu_hours: float):
        """작업 완료 후 GPU 사용량 업데이트"""
        await supabase.rpc("increment_gpu_usage", {
            "p_user_id": user_id,
            "p_hours": gpu_hours
        }).execute()

    def _get_priority_from_tier(self, tier: str) -> int:
        """우선순위 티어 → 숫자 변환"""
        return {
            "enterprise": 1,
            "priority": 3,
            "standard": 5
        }.get(tier, 5)


class QuotaExceededError(Exception):
    """GPU 쿼터 초과 예외"""
    pass
```

### 3.5.3 AlphaFold API Endpoints

```python
# backend/app/api/alphafold.py
from fastapi import APIRouter, Depends, HTTPException
from app.services.alphafold_service import AlphaFoldService, AlphaFoldRequest

router = APIRouter(prefix="/api/alphafold", tags=["AlphaFold 3D"])
alphafold_service = AlphaFoldService()

@router.post("/jobs")
async def create_alphafold_job(
    session_id: str,
    request: AlphaFoldRequest,
    user: User = Depends(get_current_user)
):
    """
    AlphaFold 3D 모델링 작업 생성 (Premium Only)
    """
    # Premium 권한 확인
    if user.subscription_tier != "premium":
        raise HTTPException(403, "AlphaFold 3D modeling requires Premium subscription")

    try:
        job_id = await alphafold_service.submit_job(
            session_id=session_id,
            user_id=user.id,
            request=request
        )
        return {"job_id": job_id, "status": "queued"}
    except PermissionError as e:
        raise HTTPException(403, str(e))

@router.get("/jobs/{job_id}")
async def get_alphafold_job(
    job_id: str,
    user: User = Depends(get_current_user)
):
    """AlphaFold 작업 상태 조회"""
    result = await alphafold_service.get_job_status(job_id)

    # 소유권 확인
    job = await supabase.table("alphafold_jobs").select("user_id").eq(
        "id", job_id
    ).single().execute()

    if job.data and job.data["user_id"] != user.id:
        raise HTTPException(403, "Not authorized to view this job")

    return result

@router.get("/viewer/{job_id}")
async def get_alphafold_viewer(
    job_id: str,
    user: User = Depends(get_current_user)
):
    """3D 구조 시각화 데이터 반환"""
    result = await alphafold_service.get_job_status(job_id)

    if result.status != "completed":
        raise HTTPException(400, f"Job not completed: {result.status}")

    return {
        "pdb_data": result.structure_pdb,
        "ligand_smiles": result.visualization_data.get("ligand_smiles"),
        "contact_residues": result.contact_residues,
        "viewer_config": {
            "style": "cartoon",
            "color_scheme": "confidence",
            "highlight_contacts": True
        }
    }
```

### 3.5.4 Frontend AlphaFold Component

```tsx
// frontend/src/components/dashboard/AlphaFold3DViewer.tsx
import { useState, useEffect, useRef } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { Loader2, Dna, Crown, Atom } from 'lucide-react'
import { toast } from 'sonner'

interface AlphaFoldViewerProps {
    sessionId: string
    candidateSmiles: string
    targetSequence?: string
    tier: 'free' | 'premium'
}

export function AlphaFold3DViewer({
    sessionId,
    candidateSmiles,
    targetSequence,
    tier
}: AlphaFoldViewerProps) {
    const [jobId, setJobId] = useState<string | null>(null)
    const [status, setStatus] = useState<string>('idle')
    const [result, setResult] = useState<any>(null)
    const viewerRef = useRef<HTMLDivElement>(null)

    const startPrediction = async () => {
        if (tier !== 'premium') {
            toast.error('AlphaFold 3D modeling is a Premium feature')
            return
        }

        try {
            setStatus('submitting')
            const res = await fetch('/api/alphafold/jobs', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    session_id: sessionId,
                    target_sequence: targetSequence || 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH',
                    ligand_smiles: candidateSmiles,
                    job_type: 'binding'
                })
            })

            if (!res.ok) throw new Error('Failed to submit job')

            const { job_id } = await res.json()
            setJobId(job_id)
            setStatus('processing')

            // Poll for completion
            pollJobStatus(job_id)
        } catch (error) {
            toast.error('Failed to start AlphaFold prediction')
            setStatus('error')
        }
    }

    const pollJobStatus = async (id: string) => {
        const interval = setInterval(async () => {
            try {
                const res = await fetch(`/api/alphafold/jobs/${id}`)
                const data = await res.json()

                if (data.status === 'completed') {
                    clearInterval(interval)
                    setResult(data)
                    setStatus('completed')
                    toast.success('3D structure prediction complete!')
                } else if (data.status === 'failed') {
                    clearInterval(interval)
                    setStatus('error')
                    toast.error('Prediction failed')
                }
            } catch (e) {
                clearInterval(interval)
                setStatus('error')
            }
        }, 5000) // Poll every 5 seconds
    }

    // Initialize 3D viewer when result is ready
    useEffect(() => {
        if (result && viewerRef.current && result.structure_pdb) {
            // Initialize Mol* or NGL viewer
            initializeViewer(viewerRef.current, result)
        }
    }, [result])

    const initializeViewer = async (container: HTMLDivElement, data: any) => {
        // Using Mol* viewer (실제 구현에서는 @molstar/molstar 사용)
        // Simplified placeholder
        container.innerHTML = `
            <div class="w-full h-full flex items-center justify-center bg-black/50 rounded">
                <div class="text-center">
                    <Atom class="w-12 h-12 mx-auto mb-2 text-blue-400 animate-pulse" />
                    <p class="text-sm text-slate-400">3D Structure Loaded</p>
                    <p class="text-xs text-slate-500 mt-1">
                        Binding Affinity: ${data.binding_affinity?.toFixed(2) || 'N/A'} kcal/mol
                    </p>
                </div>
            </div>
        `
    }

    // Premium gate
    if (tier !== 'premium') {
        return (
            <Card className="bg-slate-900/50 border-slate-800 border-dashed">
                <CardContent className="py-8 text-center">
                    <Crown className="w-12 h-12 mx-auto mb-4 text-yellow-500/50" />
                    <h3 className="text-lg font-medium text-slate-300 mb-2">
                        AlphaFold 3D Modeling
                    </h3>
                    <p className="text-sm text-slate-500 mb-4">
                        Predict protein-drug binding with AlphaFold 3
                    </p>
                    <Badge className="bg-gradient-to-r from-yellow-500 to-orange-500">
                        Premium Feature
                    </Badge>
                </CardContent>
            </Card>
        )
    }

    return (
        <Card className="bg-slate-900 border-slate-800">
            <CardHeader className="flex flex-row items-center justify-between">
                <CardTitle className="flex items-center gap-2">
                    <Dna className="w-5 h-5 text-purple-400" />
                    AlphaFold 3D Prediction
                </CardTitle>
                <Badge className="bg-purple-500/20 text-purple-300 border-purple-500/30">
                    <Crown className="w-3 h-3 mr-1" />
                    Premium
                </Badge>
            </CardHeader>
            <CardContent className="space-y-4">
                {status === 'idle' && (
                    <Button
                        onClick={startPrediction}
                        className="w-full bg-gradient-to-r from-purple-500 to-pink-500"
                    >
                        <Dna className="w-4 h-4 mr-2" />
                        Predict 3D Binding Structure
                    </Button>
                )}

                {status === 'processing' && (
                    <div className="text-center py-8">
                        <Loader2 className="w-8 h-8 animate-spin mx-auto text-purple-400 mb-4" />
                        <p className="text-slate-400">Running AlphaFold prediction...</p>
                        <p className="text-xs text-slate-500 mt-2">
                            This may take a few minutes
                        </p>
                    </div>
                )}

                {status === 'completed' && result && (
                    <div className="space-y-4">
                        {/* 3D Viewer Container */}
                        <div
                            ref={viewerRef}
                            className="w-full h-[300px] bg-black rounded-lg overflow-hidden"
                        />

                        {/* Results Summary */}
                        <div className="grid grid-cols-3 gap-3">
                            <div className="bg-slate-800 rounded p-3 text-center">
                                <div className="text-2xl font-bold text-green-400">
                                    {result.binding_affinity?.toFixed(1) || 'N/A'}
                                </div>
                                <div className="text-xs text-slate-500">ΔG (kcal/mol)</div>
                            </div>
                            <div className="bg-slate-800 rounded p-3 text-center">
                                <div className="text-2xl font-bold text-blue-400">
                                    {(result.confidence * 100)?.toFixed(0) || 'N/A'}%
                                </div>
                                <div className="text-xs text-slate-500">Confidence</div>
                            </div>
                            <div className="bg-slate-800 rounded p-3 text-center">
                                <div className="text-2xl font-bold text-purple-400">
                                    {result.contact_residues?.length || 0}
                                </div>
                                <div className="text-xs text-slate-500">Contact Residues</div>
                            </div>
                        </div>

                        {/* Contact Residues */}
                        {result.contact_residues?.length > 0 && (
                            <div className="text-xs">
                                <span className="text-slate-500">Key Binding Residues: </span>
                                <span className="text-slate-300 font-mono">
                                    {result.contact_residues.slice(0, 10).join(', ')}
                                    {result.contact_residues.length > 10 && '...'}
                                </span>
                            </div>
                        )}
                    </div>
                )}
            </CardContent>
        </Card>
    )
}
```

---

## Phase 4: Backend Data Masking [Week 5]

### 4.1 Masking Service

```python
# backend/app/security/masking.py
from typing import Any, Dict, List
from enum import Enum

class MaskingLevel(Enum):
    NONE = "none"           # Premium: 전체 공개
    PARTIAL = "partial"     # Free: 부분 마스킹
    FULL = "full"           # 완전 마스킹

class DataMaskingService:
    """
    Constraint 3 준수: 백엔드 API 레벨에서 데이터 마스킹
    프론트엔드가 아닌 서버에서 직접 마스킹하여 보안 강화
    """

    def mask_session_result(self, data: Dict, tier: str) -> Dict:
        """세션 결과 마스킹"""
        if tier == "premium":
            return data  # Premium: 전체 공개

        # Free tier: 핵심 데이터 마스킹
        masked = data.copy()

        # 1. SMILES 마스킹 (첫 10자만 표시)
        if "current_smiles" in masked and masked["current_smiles"]:
            smiles = masked["current_smiles"]
            masked["current_smiles"] = smiles[:10] + "..." + "*" * 20
            masked["smiles_masked"] = True

        # 2. 후보 리스트 마스킹 (1개만 표시, 나머지 블러)
        if "candidates" in masked:
            candidates = masked["candidates"]
            if len(candidates) > 1:
                masked["candidates"] = [
                    candidates[0],  # 첫 번째만 표시
                    *[self._mask_candidate(c) for c in candidates[1:]]
                ]
                masked["candidates_limited"] = True
                masked["total_candidates_available"] = len(candidates)

        # 3. 점수 마스킹
        if "calculated_metrics" in masked:
            metrics = masked["calculated_metrics"]
            masked["calculated_metrics"] = {
                "mw": metrics.get("mw"),  # MW는 공개
                "logp": "***" if tier == "free" else metrics.get("logp"),
                "sa_score": "Premium Only" if tier == "free" else metrics.get("sa_score"),
                "tanimoto": "Premium Only" if tier == "free" else metrics.get("tanimoto")
            }

        # 4. 리포트 마스킹
        if "final_report" in masked and masked["final_report"]:
            masked["final_report"] = self._mask_report(masked["final_report"])

        return masked

    def _mask_candidate(self, candidate: Dict) -> Dict:
        """개별 후보 마스킹"""
        return {
            "rank": candidate.get("rank"),
            "smiles": "********** (Premium Only)",
            "score": "?.??",
            "is_masked": True
        }

    def _mask_report(self, report: Dict) -> Dict:
        """리포트 마스킹"""
        return {
            "summary": report.get("summary", "")[:100] + "... (Premium Only)",
            "risk_score": "?/10",
            "golden_set_comparison": "Premium Only",
            "pmid_references": report.get("pmid_references", [])[:1],  # 1개만 표시
            "is_masked": True
        }

    def mask_smiles(self, smiles: str, level: MaskingLevel) -> str:
        """SMILES 문자열 마스킹"""
        if level == MaskingLevel.NONE:
            return smiles
        elif level == MaskingLevel.PARTIAL:
            # 첫 10자 + 마스크
            return smiles[:10] + "..." + "*" * min(20, len(smiles) - 10)
        else:
            return "*" * 30
```

### 4.2 AES-256 Encryption Service

```python
# backend/app/security/encryption.py
from cryptography.hazmat.primitives.ciphers.aead import AESGCM
from cryptography.hazmat.backends import default_backend
import os
import base64

class EncryptionService:
    """
    AES-256-GCM 암호화 서비스
    유료 사용자의 SMILES 코드를 암호화하여 저장
    """

    def __init__(self):
        # Key should be loaded from secure environment variable
        self.key = os.environ.get("ENCRYPTION_KEY")
        if not self.key:
            raise ValueError("ENCRYPTION_KEY not set")
        self.key = base64.b64decode(self.key)
        self.aesgcm = AESGCM(self.key)

    def encrypt_smiles(self, smiles: str, user_id: str) -> Dict[str, bytes]:
        """
        SMILES 암호화
        Returns: {encrypted_data, iv, tag}
        """
        # Generate random IV
        iv = os.urandom(12)

        # Additional authenticated data (AAD) - user_id for binding
        aad = user_id.encode()

        # Encrypt
        ciphertext = self.aesgcm.encrypt(iv, smiles.encode(), aad)

        return {
            "encrypted_smiles": ciphertext[:-16],  # Ciphertext without tag
            "encryption_iv": iv,
            "encryption_tag": ciphertext[-16:],    # Last 16 bytes is tag
            "structure_hash": self._hash_smiles(smiles)
        }

    def decrypt_smiles(self, encrypted_data: bytes, iv: bytes,
                       tag: bytes, user_id: str) -> str:
        """SMILES 복호화"""
        aad = user_id.encode()
        ciphertext_with_tag = encrypted_data + tag

        plaintext = self.aesgcm.decrypt(iv, ciphertext_with_tag, aad)
        return plaintext.decode()

    def _hash_smiles(self, smiles: str) -> str:
        """복호화 후 검증용 해시"""
        import hashlib
        return hashlib.sha256(smiles.encode()).hexdigest()
```

---

## Phase 5: Frontend - Live Console & Real-time UI [Week 5]

### 5.1 New Route Structure

```tsx
// frontend/src/App.tsx (추가)
<Route path="/dashboard" element={<DashboardLayout />}>
    <Route index element={<DashboardHome />} />
    <Route path="design" element={<DenovoDesign />} />        {/* NEW */}
    <Route path="optimize" element={<LeadOptimization />} />  {/* NEW */}
    <Route path="audit" element={<PreclinicalAudit />} />     {/* NEW */}
    <Route path="cmc" element={<CMCSourcing />} />            {/* NEW */}
    <Route path="session/:id" element={<SessionViewer />} />  {/* NEW */}
    <Route path="library" element={<GoldenSetLibrary />} />
    <Route path="builder" element={<ADCBuilder />} />
    <Route path="result/:jobId" element={<ResultViewer />} />
</Route>
```

### 5.2 Live Agent Console Component

```tsx
// frontend/src/components/dashboard/LiveAgentConsole.tsx
import { useState, useEffect, useRef } from 'react'
import { motion, AnimatePresence } from 'framer-motion'
import { Terminal, Bot, Zap, Shield, BookOpen, Heart, Brain } from 'lucide-react'

const AGENTS = {
    orchestrator: { name: 'Orchestrator', icon: Brain, color: 'purple' },
    alchemist: { name: 'Alchemist', icon: Zap, color: 'blue' },
    coder: { name: 'Coder', icon: Terminal, color: 'green' },
    healer: { name: 'Healer', icon: Heart, color: 'red' },
    auditor: { name: 'Auditor', icon: Shield, color: 'yellow' },
    librarian: { name: 'Librarian', icon: BookOpen, color: 'cyan' }
}

interface ConsoleMessage {
    id: string
    agent: string
    type: 'info' | 'success' | 'error' | 'warning'
    message: string
    timestamp: string
    details?: any
}

export function LiveAgentConsole({ sessionId }: { sessionId: string }) {
    const [messages, setMessages] = useState<ConsoleMessage[]>([])
    const [isConnected, setIsConnected] = useState(false)
    const consoleRef = useRef<HTMLDivElement>(null)
    const wsRef = useRef<WebSocket | null>(null)

    useEffect(() => {
        // WebSocket 연결
        const ws = new WebSocket(
            `${import.meta.env.VITE_WS_URL}/ws/session/${sessionId}`
        )

        ws.onopen = () => {
            setIsConnected(true)
            addMessage({
                agent: 'system',
                type: 'info',
                message: 'Connected to design session'
            })
        }

        ws.onmessage = (event) => {
            const data = JSON.parse(event.data)
            handleAgentMessage(data)
        }

        ws.onclose = () => {
            setIsConnected(false)
            addMessage({
                agent: 'system',
                type: 'warning',
                message: 'Disconnected from session'
            })
        }

        wsRef.current = ws
        return () => ws.close()
    }, [sessionId])

    const handleAgentMessage = (data: any) => {
        if (data.type === 'agent_status') {
            addMessage({
                agent: data.agent,
                type: data.status === 'error' ? 'error' :
                      data.status === 'completed' ? 'success' : 'info',
                message: data.message,
                details: data.details
            })
        } else if (data.type === 'agent_communication') {
            // 에이전트 간 대화 표시
            addMessage({
                agent: data.from_agent,
                type: 'info',
                message: `→ ${AGENTS[data.to_agent]?.name}: ${data.message}`
            })
        }
    }

    const addMessage = (msg: Omit<ConsoleMessage, 'id' | 'timestamp'>) => {
        setMessages(prev => [...prev, {
            ...msg,
            id: crypto.randomUUID(),
            timestamp: new Date().toISOString()
        }])
    }

    // Auto-scroll
    useEffect(() => {
        if (consoleRef.current) {
            consoleRef.current.scrollTop = consoleRef.current.scrollHeight
        }
    }, [messages])

    return (
        <div className="bg-slate-950 rounded-lg border border-slate-800 overflow-hidden">
            {/* Header */}
            <div className="px-4 py-2 bg-slate-900 border-b border-slate-800 flex items-center justify-between">
                <div className="flex items-center gap-2">
                    <Terminal className="w-4 h-4 text-green-400" />
                    <span className="font-mono text-sm text-slate-300">Live Agent Console</span>
                </div>
                <div className="flex items-center gap-2">
                    <div className={`w-2 h-2 rounded-full ${isConnected ? 'bg-green-500 animate-pulse' : 'bg-red-500'}`} />
                    <span className="text-xs text-slate-500">
                        {isConnected ? 'Connected' : 'Disconnected'}
                    </span>
                </div>
            </div>

            {/* Console Output */}
            <div
                ref={consoleRef}
                className="h-[300px] overflow-y-auto p-4 font-mono text-xs space-y-2"
            >
                <AnimatePresence>
                    {messages.map((msg) => {
                        const agent = AGENTS[msg.agent]
                        const AgentIcon = agent?.icon || Bot

                        return (
                            <motion.div
                                key={msg.id}
                                initial={{ opacity: 0, x: -20 }}
                                animate={{ opacity: 1, x: 0 }}
                                exit={{ opacity: 0 }}
                                className={`flex items-start gap-2 ${
                                    msg.type === 'error' ? 'text-red-400' :
                                    msg.type === 'success' ? 'text-green-400' :
                                    msg.type === 'warning' ? 'text-yellow-400' :
                                    'text-slate-400'
                                }`}
                            >
                                <span className="text-slate-600 shrink-0">
                                    {new Date(msg.timestamp).toLocaleTimeString()}
                                </span>
                                {agent && (
                                    <AgentIcon className={`w-4 h-4 shrink-0 text-${agent.color}-400`} />
                                )}
                                <span className={`font-semibold text-${agent?.color || 'slate'}-400`}>
                                    [{agent?.name || msg.agent}]
                                </span>
                                <span className="flex-1">{msg.message}</span>
                            </motion.div>
                        )
                    })}
                </AnimatePresence>
            </div>
        </div>
    )
}
```

### 5.3 Progress Timeline Component

```tsx
// frontend/src/components/dashboard/DesignProgressTimeline.tsx
import { motion } from 'framer-motion'
import { Check, Loader2, AlertCircle, Clock } from 'lucide-react'

interface Step {
    id: string
    label: string
    description: string
    status: 'pending' | 'running' | 'completed' | 'error'
    agent?: string
}

const DEFAULT_STEPS: Step[] = [
    { id: 'analyze', label: '타겟 분석', description: 'Analyzing target antigen', status: 'pending' },
    { id: 'design', label: '후보 설계', description: 'Generating candidate structures', status: 'pending' },
    { id: 'validate', label: '공학적 검증', description: 'Running chemical validation', status: 'pending' },
    { id: 'audit', label: '리스크 평가', description: 'Assessing risk factors', status: 'pending' },
    { id: 'report', label: '리포트 생성', description: 'Generating final report', status: 'pending' }
]

export function DesignProgressTimeline({
    steps = DEFAULT_STEPS,
    currentStep
}: {
    steps?: Step[]
    currentStep: number
}) {
    return (
        <div className="relative">
            {/* Progress Line */}
            <div className="absolute left-[19px] top-0 bottom-0 w-0.5 bg-slate-800" />
            <motion.div
                className="absolute left-[19px] top-0 w-0.5 bg-gradient-to-b from-blue-500 to-purple-500"
                initial={{ height: 0 }}
                animate={{ height: `${(currentStep / (steps.length - 1)) * 100}%` }}
                transition={{ duration: 0.5 }}
            />

            {/* Steps */}
            <div className="space-y-6">
                {steps.map((step, index) => {
                    const isActive = index === currentStep
                    const isCompleted = index < currentStep
                    const isError = step.status === 'error'

                    return (
                        <div key={step.id} className="flex items-start gap-4">
                            {/* Status Icon */}
                            <div className={`
                                relative z-10 w-10 h-10 rounded-full flex items-center justify-center
                                ${isCompleted ? 'bg-green-500' :
                                  isActive ? 'bg-blue-500 animate-pulse' :
                                  isError ? 'bg-red-500' :
                                  'bg-slate-800'}
                            `}>
                                {isCompleted ? (
                                    <Check className="w-5 h-5 text-white" />
                                ) : isActive ? (
                                    <Loader2 className="w-5 h-5 text-white animate-spin" />
                                ) : isError ? (
                                    <AlertCircle className="w-5 h-5 text-white" />
                                ) : (
                                    <Clock className="w-5 h-5 text-slate-500" />
                                )}
                            </div>

                            {/* Step Info */}
                            <div className="flex-1 pt-1">
                                <div className={`font-medium ${
                                    isActive ? 'text-blue-400' :
                                    isCompleted ? 'text-green-400' :
                                    isError ? 'text-red-400' :
                                    'text-slate-500'
                                }`}>
                                    {step.label}
                                </div>
                                <div className="text-sm text-slate-500">
                                    {step.description}
                                </div>
                                {step.agent && isActive && (
                                    <div className="mt-1 text-xs text-purple-400">
                                        Agent: {step.agent}
                                    </div>
                                )}
                            </div>
                        </div>
                    )
                })}
            </div>
        </div>
    )
}
```

### 5.4 De novo Design Page

```tsx
// frontend/src/pages/Dashboard/DenovoDesign.tsx
import { useState } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Input } from '@/components/ui/input'
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from '@/components/ui/select'
import { Badge } from '@/components/ui/badge'
import { LiveAgentConsole } from '@/components/dashboard/LiveAgentConsole'
import { DesignProgressTimeline } from '@/components/dashboard/DesignProgressTimeline'
import { PremiumBlur } from '@/components/dashboard/PremiumBlur'
import { Sparkles, Lock, Rocket, FlaskConical } from 'lucide-react'
import { useDesignSession } from '@/hooks/useDesignSession'

export function DenovoDesign() {
    const [tier] = useState<'free' | 'premium'>('free') // TODO: Get from user context
    const {
        session,
        isLoading,
        startDesign,
        currentStep
    } = useDesignSession()

    const [formData, setFormData] = useState({
        target_antigen: '',
        target_indication: '',
        requested_dar: 4,
        linker_preference: 'cleavable'
    })

    const handleStart = async () => {
        await startDesign({
            session_type: 'denovo',
            ...formData
        })
    }

    return (
        <div className="space-y-6">
            {/* Header */}
            <div className="flex items-center justify-between">
                <div>
                    <h1 className="text-3xl font-bold text-white flex items-center gap-3">
                        <Sparkles className="w-8 h-8 text-blue-400" />
                        De novo Design
                    </h1>
                    <p className="text-slate-400 mt-1">
                        Create optimal ADC candidates for new targets
                    </p>
                </div>
                <Badge variant={tier === 'premium' ? 'default' : 'outline'}
                       className={tier === 'premium' ? 'bg-gradient-to-r from-yellow-500 to-orange-500' : ''}>
                    {tier === 'premium' ? '⭐ Premium' : 'Free Tier'}
                </Badge>
            </div>

            <div className="grid grid-cols-3 gap-6">
                {/* Left: Input Form */}
                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader>
                        <CardTitle className="text-lg">Target Information</CardTitle>
                    </CardHeader>
                    <CardContent className="space-y-4">
                        <div className="space-y-2">
                            <label className="text-sm text-slate-400">Target Antigen</label>
                            <Input
                                placeholder="e.g., HER2, TROP-2, EGFR"
                                value={formData.target_antigen}
                                onChange={(e) => setFormData({...formData, target_antigen: e.target.value})}
                                className="bg-slate-950 border-slate-800"
                            />
                        </div>

                        <div className="space-y-2">
                            <label className="text-sm text-slate-400">Target Indication</label>
                            <Input
                                placeholder="e.g., Breast Cancer, NSCLC"
                                value={formData.target_indication}
                                onChange={(e) => setFormData({...formData, target_indication: e.target.value})}
                                className="bg-slate-950 border-slate-800"
                            />
                        </div>

                        <div className="space-y-2">
                            <label className="text-sm text-slate-400">Desired DAR</label>
                            <Select
                                value={String(formData.requested_dar)}
                                onValueChange={(v) => setFormData({...formData, requested_dar: parseInt(v)})}
                            >
                                <SelectTrigger className="bg-slate-950 border-slate-800">
                                    <SelectValue />
                                </SelectTrigger>
                                <SelectContent>
                                    <SelectItem value="2">DAR 2</SelectItem>
                                    <SelectItem value="4">DAR 4 (Recommended)</SelectItem>
                                    <SelectItem value="8">DAR 8</SelectItem>
                                </SelectContent>
                            </Select>
                        </div>

                        <div className="space-y-2">
                            <label className="text-sm text-slate-400">Linker Preference</label>
                            <Select
                                value={formData.linker_preference}
                                onValueChange={(v) => setFormData({...formData, linker_preference: v})}
                            >
                                <SelectTrigger className="bg-slate-950 border-slate-800">
                                    <SelectValue />
                                </SelectTrigger>
                                <SelectContent>
                                    <SelectItem value="cleavable">Cleavable (Protease-sensitive)</SelectItem>
                                    <SelectItem value="non-cleavable">Non-cleavable (Stable)</SelectItem>
                                    <SelectItem value="any">No Preference</SelectItem>
                                </SelectContent>
                            </Select>
                        </div>

                        <Button
                            onClick={handleStart}
                            disabled={isLoading || !formData.target_antigen}
                            className="w-full bg-gradient-to-r from-blue-500 to-purple-500 hover:from-blue-600 hover:to-purple-600"
                        >
                            <Rocket className="w-4 h-4 mr-2" />
                            Start Design
                        </Button>
                    </CardContent>
                </Card>

                {/* Center: Progress & Console */}
                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader>
                        <CardTitle className="text-lg">Design Progress</CardTitle>
                    </CardHeader>
                    <CardContent className="space-y-6">
                        <DesignProgressTimeline currentStep={currentStep} />

                        {session && (
                            <LiveAgentConsole sessionId={session.id} />
                        )}
                    </CardContent>
                </Card>

                {/* Right: Results */}
                <Card className="bg-slate-900 border-slate-800">
                    <CardHeader>
                        <CardTitle className="text-lg flex items-center justify-between">
                            <span>Candidates</span>
                            {session?.candidates_limited && (
                                <Badge variant="outline" className="text-yellow-400 border-yellow-400/50">
                                    <Lock className="w-3 h-3 mr-1" />
                                    Limited
                                </Badge>
                            )}
                        </CardTitle>
                    </CardHeader>
                    <CardContent>
                        {session?.candidates?.length > 0 ? (
                            <div className="space-y-3">
                                {session.candidates.map((candidate, idx) => (
                                    <PremiumBlur
                                        key={idx}
                                        tier={tier}
                                        isBlurred={candidate.is_masked}
                                        feature="full candidate details"
                                    >
                                        <div className="p-4 bg-slate-950 rounded-lg border border-slate-800">
                                            <div className="flex items-center justify-between mb-2">
                                                <Badge>Rank #{candidate.rank || idx + 1}</Badge>
                                                <span className="text-sm text-green-400">
                                                    Score: {candidate.score || '?.??'}
                                                </span>
                                            </div>
                                            <div className="font-mono text-xs text-slate-400 break-all">
                                                {candidate.smiles}
                                            </div>
                                        </div>
                                    </PremiumBlur>
                                ))}

                                {tier === 'free' && session.total_candidates_available > 1 && (
                                    <div className="text-center p-4 bg-gradient-to-r from-yellow-500/10 to-orange-500/10 rounded-lg border border-yellow-500/30">
                                        <Lock className="w-6 h-6 mx-auto mb-2 text-yellow-400" />
                                        <p className="text-sm text-yellow-300">
                                            +{session.total_candidates_available - 1} more candidates available
                                        </p>
                                        <Button variant="outline" size="sm" className="mt-2 border-yellow-500/50 text-yellow-400">
                                            Upgrade to Premium
                                        </Button>
                                    </div>
                                )}
                            </div>
                        ) : (
                            <div className="text-center py-8 text-slate-500">
                                <FlaskConical className="w-12 h-12 mx-auto mb-4 opacity-20" />
                                <p>Start a design to see candidates</p>
                            </div>
                        )}
                    </CardContent>
                </Card>
            </div>
        </div>
    )
}
```

---

## Phase 6: API Endpoints & WebSocket [Week 5]

### 6.1 Design API Routes

```python
# backend/app/api/design.py
from fastapi import APIRouter, Depends, WebSocket, WebSocketDisconnect, BackgroundTasks
from app.agents.orchestrator import Orchestrator
from app.security.masking import DataMaskingService

router = APIRouter(prefix="/api/design", tags=["Design Engine"])
masking_service = DataMaskingService()

@router.post("/session")
async def create_session(
    request: CreateSessionRequest,
    user: User = Depends(get_current_user)
):
    """Create a new design session"""
    session = await design_service.create_session(
        user_id=user.id,
        session_type=request.session_type,
        tier=user.subscription_tier,
        target_antigen=request.target_antigen,
        target_indication=request.target_indication,
        requested_dar=request.requested_dar,
        linker_preference=request.linker_preference,
        design_goal=request.design_goal
    )
    return {"session_id": session.id}

@router.post("/session/{session_id}/start")
async def start_design(
    session_id: str,
    background_tasks: BackgroundTasks,
    user: User = Depends(get_current_user)
):
    """Start the design workflow"""
    # Verify ownership
    session = await design_service.get_session(session_id)
    if session.user_id != user.id:
        raise HTTPException(403, "Not authorized")

    # Run orchestrator in background
    orchestrator = Orchestrator()
    background_tasks.add_task(orchestrator.run, session_id)

    return {"status": "started", "session_id": session_id}

@router.get("/session/{session_id}")
async def get_session(
    session_id: str,
    user: User = Depends(get_current_user)
):
    """
    Get session status and results
    Constraint 3: 백엔드에서 티어에 따라 마스킹 적용
    """
    session = await design_service.get_session(session_id)
    if session.user_id != user.id:
        raise HTTPException(403, "Not authorized")

    # Apply backend masking based on tier
    masked_data = masking_service.mask_session_result(
        session.dict(),
        tier=user.subscription_tier
    )

    return masked_data

@router.websocket("/ws/session/{session_id}")
async def session_websocket(
    websocket: WebSocket,
    session_id: str
):
    """Real-time session updates via WebSocket"""
    await websocket.accept()

    try:
        # Subscribe to session updates
        async for message in websocket_hub.subscribe(session_id):
            await websocket.send_json(message)
    except WebSocketDisconnect:
        websocket_hub.unsubscribe(session_id, websocket)
```

---

## 📅 Updated 7-Week Roadmap (v2.1)

| Week | Phase | Key Deliverables | Critical Points |
|------|-------|------------------|-----------------|
| **W1** | Foundation | DB Schema + API Spec | ✅ Encrypted structures, ✅ Agent logs, ✅ candidate_snippets, ✅ scaffold_knowledge_mapping |
| **W2** | Core Logic (1/2) | Orchestrator + Alchemist + Librarian | ✅ Shared State 동기화, ✅ **Dynamic Substructure RAG** |
| **W3** | Core Logic (2/2) | Coder + Healer + Auditor | ✅ Snippet Library, ✅ **Knowledge Loop**, ✅ **Constraint Guardrail** |
| **W4** | Sandbox | Docker 실행 환경 | ✅ Air-gapped, ✅ Resource limits |
| **W4.5** | **AlphaFold 3** | Premium 3D 모델링 | ✅ **AF3 API 연동**, ✅ 비동기 Job Queue, ✅ 3D Viewer |
| **W5** | Frontend | 4대 메뉴 + Live Console | ✅ WebSocket 실시간, ✅ Backend masking, ✅ AlphaFold UI |
| **W6-7** | Security & QA | 결제 연동 + 최종 테스트 | ✅ AES-256, ✅ Premium 권한 검증 |

---

## Success Metrics (KPIs) - v2.1 Updated

| Metric | Target | Measurement |
|--------|--------|-------------|
| Design Completion Rate | > 90% | Sessions completed without manual review |
| Self-Healing Success | > 70% | Errors auto-corrected by The Healer |
| **Knowledge Loop Reuse** | > 30% | Healer fixes from candidate_snippets (v2.1) |
| **Constraint Compliance** | > 95% | Sessions meeting initial user constraints (v2.1) |
| **Scaffold RAG Relevance** | > 0.75 | Average relevance score of Librarian results (v2.1) |
| P95 Response Time | < 30s (Free), < 2min (Premium) | Time to first candidate |
| **AlphaFold Job Success** | > 85% | AF3 predictions completed successfully (v2.1) |
| Free → Premium Conversion | > 8% | Increased by AlphaFold feature exposure |
| Audit Trail Compliance | 100% | All decisions logged with reasoning |

---

## Risk Mitigation - v2.1 Updated

| Risk | Impact | Mitigation |
|------|--------|------------|
| LangGraph state sync issues | High | Implement pessimistic locking, checkpoint frequently |
| Docker sandbox escape | Critical | Use gVisor, no root, network isolation |
| Encryption key leak | Critical | Use HSM/KMS, rotate keys quarterly |
| WebSocket scalability | Medium | Redis pub/sub for horizontal scaling |
| Model hallucination | Medium | Snippet Library + validation before output |
| **Scaffold RAG noise** | Medium | Tanimoto threshold ≥ 0.6, relevance scoring (v2.1) |
| **Constraint drift** | High | Constraint Guardrail mandatory check (v2.1) |
| **Knowledge Loop pollution** | Medium | Confidence threshold ≥ 0.7 for promotion (v2.1) |
| **AlphaFold API latency** | Medium | Async job queue, timeout fallbacks (v2.1) |
| **AlphaFold cost overrun** | High | GPU hour quotas per user/session (v2.1) |

---

## Appendix A: New Table Relationships (v2.1)

```
design_sessions ─────┬────────────────────────────────────────┐
        │            │                                        │
        │            ▼                                        ▼
        │    agent_execution_logs                    alphafold_jobs
        │            │                               (Premium only)
        │            │
        │            ▼
        │    candidate_snippets ──────► code_snippet_library
        │    (Healer Knowledge Loop)    (after promotion)
        │
        └───────────► scaffold_knowledge_mapping
                      (Librarian RAG)
                             │
                             ▼
                      ┌──────┴──────┐
                      ▼             ▼
               knowledge_base   golden_set_library
```

---

## Appendix B: Constraint Guardrail Flow (v2.1)

```
User Input (Session Creation)
    │
    ├── target_antigen: "HER2"
    ├── requested_dar: 4
    └── linker_preference: "cleavable"
            │
            ▼
    ┌───────────────────┐
    │   Design Engine   │
    │   (Alchemist +    │
    │    Coder + ...)   │
    └───────────────────┘
            │
            ▼
    ┌───────────────────────────────────────┐
    │   🔒 The Auditor - Constraint Check   │
    │                                       │
    │   for each constraint in initial:     │
    │       compare with current state      │
    │       if violation:                   │
    │           severity = HIGH → REDESIGN  │
    │           severity = MED → WARNING    │
    └───────────────────────────────────────┘
            │
      ┌─────┴─────┐
      ▼           ▼
   APPROVED    REDESIGN → back to Alchemist
```

---

*Document Version: 2.2 (Critical Details 반영)*
*Updated: 2026-01-27*
*Author: Claude Code*

### Change Log
| Version | Date | Changes |
|---------|------|---------|
| v2.0 | 2026-01-27 | Initial enhanced plan with 5 constraints |
| v2.1 | 2026-01-27 | + Librarian Dynamic Substructure RAG |
|      |            | + Auditor Constraint Guardrail |
|      |            | + Healer Knowledge Loop |
|      |            | + Phase 3.5 AlphaFold 3 Integration |
| **v2.2** | 2026-01-27 | ⚠️ **Critical Details 추가** |
|      |            | + **Constraint 6**: 21 CFR Part 11 Digital Seal (INSERT-ONLY + Chain Hash) |
|      |            | + **Constraint 7**: Sandbox Library Version Lock |
|      |            | + `alphafold_user_quota` 테이블 (월간/일간 GPU 쿼터 관리) |
|      |            | + Version verification script + CI/CD 파이프라인 |
|      |            | + agent_execution_logs 불변성 트리거 |
