-- ADC-GenAI Platform Database Schema
-- Supabase (PostgreSQL) 스키마
-- Version: 2.1 (RLS Policy Fixes)

-- ============================================================
-- 익스텐션 설정 (Supabase Dashboard에서 먼저 설정 필요!)
-- Database > Extensions > vector 검색 후 Enable
-- ============================================================
CREATE EXTENSION IF NOT EXISTS vector;


-- ============================================================
-- 사용자 프로필 테이블 (Supabase Auth와 연동)
-- ============================================================
CREATE TABLE IF NOT EXISTS profiles (
    id UUID PRIMARY KEY REFERENCES auth.users(id) ON DELETE CASCADE,
    email TEXT UNIQUE NOT NULL,
    full_name TEXT,
    organization TEXT,
    role TEXT DEFAULT 'user' CHECK (role IN ('user', 'admin', 'super_admin')),
    plan TEXT DEFAULT 'free' CHECK (plan IN ('free', 'pro', 'enterprise')),
    credits INTEGER DEFAULT 50,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- ============================================================
-- 크레딧 트랜잭션 테이블 (매출 정산용)
-- ============================================================
CREATE TABLE IF NOT EXISTS credit_transactions (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    user_id UUID REFERENCES profiles(id) ON DELETE CASCADE,
    amount INTEGER NOT NULL, -- 양수: 지급, 음수: 차감
    balance_after INTEGER NOT NULL, -- 거래 후 잔액
    transaction_type TEXT NOT NULL CHECK (transaction_type IN (
        'grant',          -- 관리자가 수동 지급
        'purchase',       -- 결제로 구매
        'simulation',     -- 시뮬레이션 사용
        'refund',         -- 환불
        'expiry'          -- 만료
    )),
    reason TEXT,          -- 지급 사유 (영업 선물 등)
    granted_by UUID REFERENCES profiles(id), -- 관리자 ID (grant인 경우)
    job_id TEXT,          -- 관련 작업 ID (simulation인 경우)
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- ============================================================
-- 시뮬레이션 작업 테이블
-- ============================================================
CREATE TABLE IF NOT EXISTS jobs (
    id TEXT PRIMARY KEY,
    user_id UUID REFERENCES profiles(id) ON DELETE CASCADE,
    status TEXT DEFAULT 'queued' CHECK (status IN (
        'queued', 'running', 'completed', 'failed', 'partial'
    )),
    progress INTEGER DEFAULT 0,
    mode TEXT DEFAULT 'deep' CHECK (mode IN ('fast', 'deep')),
    credits_used INTEGER DEFAULT 0,
    
    -- 입력 데이터
    input_data JSONB NOT NULL,
    
    -- 에이전트 상태
    agent_statuses JSONB DEFAULT '{
        "structure": "pending",
        "toxicology": "pending",
        "patent": "pending",
        "competitor": "pending",
        "clinical": "pending",
        "report": "pending"
    }'::jsonb,
    
    -- 분석 결과
    structure_analysis JSONB,
    toxicity_risks JSONB,
    patent_landscape JSONB,
    competitors JSONB,
    clinical_protocol JSONB,
    
    -- 최종 결과
    final_grade TEXT,
    recommendation TEXT,
    executive_summary TEXT,
    scores JSONB,
    report_url TEXT,
    
    -- 에러 로그
    errors JSONB DEFAULT '[]'::jsonb,
    
    -- 타임스탬프
    created_at TIMESTAMPTZ DEFAULT NOW(),
    started_at TIMESTAMPTZ,
    completed_at TIMESTAMPTZ
);

-- ============================================================
-- Golden Set 테이블 (FDA 승인 ADC 데이터)
-- [Human-in-the-Loop] 수집 -> 검토 -> 배포 워크플로우
-- ============================================================
CREATE TABLE IF NOT EXISTS golden_set (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    
    -- 기본 정보
    drug_name TEXT NOT NULL,
    brand_name TEXT,
    target TEXT NOT NULL,
    antibody TEXT NOT NULL,
    payload TEXT NOT NULL,
    linker TEXT,
    dar NUMERIC(3,1),
    company TEXT,
    approval_date DATE,
    indications TEXT[],
    
    -- 상세 데이터
    clinical_data JSONB,
    structure_data JSONB,
    
    -- =====================================================
    -- [Human-in-the-Loop] 상태 관리 컬럼
    -- =====================================================
    status TEXT DEFAULT 'draft' CHECK (status IN ('draft', 'approved', 'rejected')),
    enrichment_source TEXT,  -- e.g., 'perplexity_sonar_medium', 'clinical_trials', 'pubmed'
    reviewer_id UUID REFERENCES profiles(id), -- 승인/반려한 관리자
    reviewer_note TEXT,      -- 관리자 반려/수정 사유
    raw_data JSONB,          -- 원본 데이터 보관용 (디버깅)
    approved_at TIMESTAMPTZ, -- 승인 시간
    
    -- 타임스탬프
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Golden Set 벡터 임베딩 테이블 (RAG용)
CREATE TABLE IF NOT EXISTS golden_set_embeddings (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    golden_set_id UUID REFERENCES golden_set(id) ON DELETE CASCADE,
    content TEXT NOT NULL,  -- 임베딩 생성용 텍스트
    embedding VECTOR(1536), -- OpenAI ada-002 임베딩
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- ============================================================
-- AI 에이전트 프롬프트 테이블 (버전 관리)
-- ============================================================
CREATE TABLE IF NOT EXISTS agent_prompts (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    agent_id TEXT NOT NULL, -- structure, toxicology, patent, etc.
    version TEXT NOT NULL,
    is_live BOOLEAN DEFAULT FALSE,
    system_prompt TEXT NOT NULL,
    created_by UUID REFERENCES profiles(id),
    created_at TIMESTAMPTZ DEFAULT NOW(),
    published_at TIMESTAMPTZ,
    
    UNIQUE(agent_id, version)
);

-- ============================================================
-- 데이터 소스 동기화 로그
-- ============================================================
CREATE TABLE IF NOT EXISTS data_sync_logs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    source_id TEXT NOT NULL, -- clinical, pubmed, news, goldenset
    status TEXT DEFAULT 'running' CHECK (status IN ('running', 'completed', 'failed')),
    records_synced INTEGER DEFAULT 0,
    records_drafted INTEGER DEFAULT 0,  -- draft 상태로 저장된 수
    error_message TEXT,
    started_at TIMESTAMPTZ DEFAULT NOW(),
    completed_at TIMESTAMPTZ
);

-- ============================================================
-- 인덱스
-- ============================================================
CREATE INDEX IF NOT EXISTS idx_jobs_user_id ON jobs(user_id);
CREATE INDEX IF NOT EXISTS idx_jobs_status ON jobs(status);
CREATE INDEX IF NOT EXISTS idx_jobs_created_at ON jobs(created_at DESC);
CREATE INDEX IF NOT EXISTS idx_transactions_user_id ON credit_transactions(user_id);
CREATE INDEX IF NOT EXISTS idx_transactions_created_at ON credit_transactions(created_at DESC);
CREATE INDEX IF NOT EXISTS idx_golden_set_target ON golden_set(target);

-- [Human-in-the-Loop] 관리자 페이지 조회 속도 향상
CREATE INDEX IF NOT EXISTS idx_golden_set_status ON golden_set(status);
CREATE INDEX IF NOT EXISTS idx_golden_set_created_at ON golden_set(created_at DESC);

-- ============================================================
-- RLS (Row Level Security) 정책
-- ============================================================
ALTER TABLE profiles ENABLE ROW LEVEL SECURITY;
ALTER TABLE jobs ENABLE ROW LEVEL SECURITY;
ALTER TABLE credit_transactions ENABLE ROW LEVEL SECURITY;
ALTER TABLE golden_set ENABLE ROW LEVEL SECURITY;
ALTER TABLE golden_set_embeddings ENABLE ROW LEVEL SECURITY;

-- 사용자는 자신의 프로필만 조회/수정 가능
DROP POLICY IF EXISTS "Users can view own profile" ON profiles;
CREATE POLICY "Users can view own profile" ON profiles
    FOR SELECT USING (auth.uid() = id);

DROP POLICY IF EXISTS "Users can update own profile" ON profiles;
CREATE POLICY "Users can update own profile" ON profiles
    FOR UPDATE USING (auth.uid() = id);

-- 사용자는 자신의 작업만 조회 가능
DROP POLICY IF EXISTS "Users can view own jobs" ON jobs;
CREATE POLICY "Users can view own jobs" ON jobs
    FOR SELECT USING (auth.uid() = user_id);

DROP POLICY IF EXISTS "Users can create own jobs" ON jobs;
CREATE POLICY "Users can create own jobs" ON jobs
    FOR INSERT WITH CHECK (auth.uid() = user_id);

-- =====================================================
-- [Human-in-the-Loop] Golden Set RLS 정책
-- AI Agents (RAG Reader): 오직 'approved' 상태만 조회 가능
-- Admin: 모든 상태 조회 가능
-- =====================================================
DROP POLICY IF EXISTS "Agents read only approved data" ON golden_set;
CREATE POLICY "Agents read only approved data" ON golden_set
    FOR SELECT USING (
        status = 'approved' 
        OR auth.role() = 'service_role'
        OR EXISTS (SELECT 1 FROM profiles WHERE id = auth.uid() AND role IN ('admin', 'super_admin'))
    );

DROP POLICY IF EXISTS "Admins can manage golden_set" ON golden_set;
CREATE POLICY "Admins can manage golden_set" ON golden_set
    FOR ALL USING (
        EXISTS (SELECT 1 FROM profiles WHERE id = auth.uid() AND role IN ('admin', 'super_admin'))
        OR auth.role() = 'service_role'
    );

-- 임베딩도 approved 데이터만 조회 가능
DROP POLICY IF EXISTS "Read approved embeddings" ON golden_set_embeddings;
CREATE POLICY "Read approved embeddings" ON golden_set_embeddings
    FOR SELECT USING (
        EXISTS (
            SELECT 1 FROM golden_set 
            WHERE golden_set.id = golden_set_embeddings.golden_set_id 
            AND golden_set.status = 'approved'
        )
        OR auth.role() = 'service_role'
    );

-- 관리자는 모든 데이터 접근 가능
DROP POLICY IF EXISTS "Admins can access all profiles" ON profiles;
CREATE POLICY "Admins can access all profiles" ON profiles
    FOR ALL USING (
        EXISTS (SELECT 1 FROM profiles WHERE id = auth.uid() AND role IN ('admin', 'super_admin'))
    );

DROP POLICY IF EXISTS "Admins can access all jobs" ON jobs;
CREATE POLICY "Admins can access all jobs" ON jobs
    FOR ALL USING (
        EXISTS (SELECT 1 FROM profiles WHERE id = auth.uid() AND role IN ('admin', 'super_admin'))
    );

-- ============================================================
-- 트리거: updated_at 자동 갱신
-- ============================================================
CREATE OR REPLACE FUNCTION update_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS profiles_updated_at ON profiles;
CREATE TRIGGER profiles_updated_at
    BEFORE UPDATE ON profiles
    FOR EACH ROW EXECUTE FUNCTION update_updated_at();

DROP TRIGGER IF EXISTS golden_set_updated_at ON golden_set;
CREATE TRIGGER golden_set_updated_at
    BEFORE UPDATE ON golden_set
    FOR EACH ROW EXECUTE FUNCTION update_updated_at();

-- ============================================================
-- [Human-in-the-Loop] 승인 시 자동 임베딩 생성 함수
-- 이 함수는 Edge Function이나 API에서 호출
-- ============================================================
CREATE OR REPLACE FUNCTION approve_golden_set(
    p_golden_set_id UUID,
    p_reviewer_id UUID,
    p_note TEXT DEFAULT NULL
)
RETURNS VOID AS $$
BEGIN
    -- 1. 상태를 approved로 변경
    UPDATE golden_set 
    SET 
        status = 'approved',
        reviewer_id = p_reviewer_id,
        reviewer_note = p_note,
        approved_at = NOW()
    WHERE id = p_golden_set_id;
    
    -- 2. 임베딩 생성은 API에서 OpenAI 호출 후 별도 INSERT
    -- (Supabase Edge Function 또는 Backend API에서 처리)
END;
$$ LANGUAGE plpgsql;

-- ============================================================
-- [RAG] 벡터 유사도 검색 함수 (Cosine Similarity)
-- ============================================================
CREATE OR REPLACE FUNCTION match_golden_set (
  query_embedding vector(1536),
  match_threshold float,
  match_count int
)
RETURNS TABLE (
  id uuid,
  content text,
  similarity float
)
LANGUAGE plpgsql
AS $$
BEGIN
  RETURN QUERY
  SELECT
    golden_set_embeddings.source_id,
    golden_set_embeddings.chunk_content,
    1 - (golden_set_embeddings.embedding <=> query_embedding) AS similarity
  FROM golden_set_embeddings
  WHERE 1 - (golden_set_embeddings.embedding <=> query_embedding) > match_threshold
  ORDER BY similarity DESC
  LIMIT match_count;
END;
$$;
