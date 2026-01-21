-- Job Locks 테이블 (중복 작업 방지용)
CREATE TABLE IF NOT EXISTS job_locks (
    id uuid NOT NULL DEFAULT gen_random_uuid(),
    job_type text NOT NULL UNIQUE,
    acquired_at timestamp with time zone NOT NULL DEFAULT now(),
    CONSTRAINT job_locks_pkey PRIMARY KEY (id)
);

-- 인덱스
CREATE INDEX IF NOT EXISTS idx_job_locks_type ON job_locks(job_type);
CREATE INDEX IF NOT EXISTS idx_job_locks_acquired ON job_locks(acquired_at);

-- LLM 비용 추적 테이블
CREATE TABLE IF NOT EXISTS llm_usage_logs (
    id uuid NOT NULL DEFAULT gen_random_uuid(),
    model text NOT NULL,
    input_tokens integer NOT NULL DEFAULT 0,
    output_tokens integer NOT NULL DEFAULT 0,
    cost_usd numeric(10, 6) NOT NULL DEFAULT 0,
    created_at timestamp with time zone NOT NULL DEFAULT now(),
    CONSTRAINT llm_usage_logs_pkey PRIMARY KEY (id)
);

-- 일별 비용 조회를 위한 인덱스
CREATE INDEX IF NOT EXISTS idx_llm_usage_created ON llm_usage_logs(created_at);

-- System Config 테이블 (AI Refiner 시스템 상태 관리)
CREATE TABLE IF NOT EXISTS system_config (
    key text PRIMARY KEY,
    value text NOT NULL,
    updated_at timestamp with time zone DEFAULT now()
);

-- 기본 상태 삽입
INSERT INTO system_config (key, value) VALUES ('AI_REFINER_STATUS', 'ACTIVE')
ON CONFLICT (key) DO NOTHING;

