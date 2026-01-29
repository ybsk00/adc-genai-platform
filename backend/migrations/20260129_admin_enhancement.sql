-- ============================================================================
-- Migration: Admin Dashboard Enhancement v2.2+
-- Date: 2026-01-29
-- Description: 관리자 대시보드 고도화 - 예산 관리, 감사 로그, AI 정제
-- ============================================================================

-- ============================================================================
-- 1. System Configuration Table (기존 테이블 확장)
-- ============================================================================

-- 기존 테이블에 누락된 컬럼 추가 (이미 있으면 무시)
DO $$
BEGIN
    -- description 컬럼 추가
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'system_config' AND column_name = 'description') THEN
        ALTER TABLE system_config ADD COLUMN description TEXT;
    END IF;

    -- updated_by 컬럼 추가
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'system_config' AND column_name = 'updated_by') THEN
        ALTER TABLE system_config ADD COLUMN updated_by UUID;
    END IF;
EXCEPTION
    WHEN undefined_table THEN
        -- 테이블이 없으면 새로 생성
        CREATE TABLE system_config (
            key TEXT PRIMARY KEY,
            value JSONB NOT NULL,
            description TEXT,
            updated_at TIMESTAMPTZ DEFAULT NOW(),
            updated_by UUID
        );
END $$;

-- 기본 설정값 (description 없이 삽입)
INSERT INTO system_config (key, value) VALUES
    ('USE_NIM_API', '"true"'),
    ('FALLBACK_ENABLED', '"true"'),
    ('EMERGENCY_STOP', '"false"')
ON CONFLICT (key) DO NOTHING;


-- ============================================================================
-- 2. Budget Configuration Table (API 비용 관리)
-- ============================================================================

CREATE TABLE IF NOT EXISTS budget_configs (
    api_provider TEXT PRIMARY KEY,
    daily_limit_usd NUMERIC(10,2) NOT NULL,
    current_usage_usd NUMERIC(10,2) DEFAULT 0,
    monthly_limit_usd NUMERIC(10,2),
    current_monthly_usd NUMERIC(10,2) DEFAULT 0,
    auto_fallback_enabled BOOLEAN DEFAULT TRUE,
    emergency_stop_enabled BOOLEAN DEFAULT FALSE,
    alert_threshold_pct INT DEFAULT 80,
    alert_email TEXT,
    alert_slack_webhook TEXT,
    last_alert_sent_at TIMESTAMPTZ,
    last_daily_reset_at TIMESTAMPTZ DEFAULT NOW(),
    last_monthly_reset_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- 초기 예산 설정
INSERT INTO budget_configs (api_provider, daily_limit_usd, monthly_limit_usd) VALUES
    ('nvidia_nim', 500.00, 10000.00),
    ('gemini', 200.00, 4000.00),
    ('alphafold', 300.00, 6000.00)
ON CONFLICT (api_provider) DO NOTHING;


-- ============================================================================
-- 3. API Usage Logs Table (API 호출 로그)
-- ============================================================================

CREATE TABLE IF NOT EXISTS api_usage_logs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    api_provider TEXT NOT NULL,
    endpoint TEXT,
    model_name TEXT,
    tokens_input INT,
    tokens_output INT,
    cost_usd NUMERIC(10,6),
    latency_ms INT,
    success BOOLEAN DEFAULT TRUE,
    error_message TEXT,
    session_id UUID,
    user_id UUID,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_api_usage_provider ON api_usage_logs(api_provider);
CREATE INDEX IF NOT EXISTS idx_api_usage_created ON api_usage_logs(created_at DESC);
CREATE INDEX IF NOT EXISTS idx_api_usage_session ON api_usage_logs(session_id);


-- ============================================================================
-- 4. Admin Audit Logs Table (관리자 감사 로그)
-- ============================================================================

CREATE TABLE IF NOT EXISTS admin_audit_logs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    admin_id UUID,
    admin_email TEXT,
    action_type TEXT NOT NULL,
    action_category TEXT,
    target_type TEXT,
    target_id TEXT,
    before_value JSONB,
    after_value JSONB,
    reason TEXT,
    ip_address TEXT,
    user_agent TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_audit_logs_admin ON admin_audit_logs(admin_id);
CREATE INDEX IF NOT EXISTS idx_audit_logs_action ON admin_audit_logs(action_type);
CREATE INDEX IF NOT EXISTS idx_audit_logs_created ON admin_audit_logs(created_at DESC);
CREATE INDEX IF NOT EXISTS idx_audit_logs_category ON admin_audit_logs(action_category);


-- ============================================================================
-- 5. Agent Prompts Version Control (SKIP - 기존 테이블 사용)
-- ============================================================================
-- 기존 agent_prompts 테이블이 agent_id 컬럼을 사용하므로 스킵
-- 필요시 별도로 ALTER TABLE로 컬럼 추가


-- ============================================================================
-- 6. Target Synonyms Table (타겟 동의어 사전)
-- ============================================================================

-- 테이블이 없으면 생성
CREATE TABLE IF NOT EXISTS target_synonyms (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    canonical_name TEXT NOT NULL UNIQUE,
    synonyms TEXT[] NOT NULL,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- 기존 테이블에 category 컬럼 추가 (없으면)
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'target_synonyms' AND column_name = 'category') THEN
        ALTER TABLE target_synonyms ADD COLUMN category TEXT;
    END IF;
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'target_synonyms' AND column_name = 'description') THEN
        ALTER TABLE target_synonyms ADD COLUMN description TEXT;
    END IF;
END $$;

-- 기본 동의어 데이터 (category 없이 삽입)
INSERT INTO target_synonyms (canonical_name, synonyms) VALUES
    ('HER2', ARRAY['ERBB2', 'neu', 'CD340', 'HER-2', 'HER2/neu']),
    ('EGFR', ARRAY['HER1', 'ERBB1', 'ErbB-1']),
    ('PD-1', ARRAY['PD1', 'PDCD1', 'CD279']),
    ('PD-L1', ARRAY['PDL1', 'CD274', 'B7-H1']),
    ('TROP2', ARRAY['TROP-2', 'TACSTD2', 'EGP-1']),
    ('CD19', ARRAY['B4', 'CVID3']),
    ('CD33', ARRAY['SIGLEC3', 'gp67']),
    ('BCMA', ARRAY['TNFRSF17', 'CD269']),
    ('VEGF', ARRAY['VEGFA', 'VPF']),
    ('CD3', ARRAY['T3'])
ON CONFLICT (canonical_name) DO NOTHING;


-- ============================================================================
-- 7. Quarantined Data Table (격리 데이터)
-- ============================================================================

CREATE TABLE IF NOT EXISTS quarantined_data (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    source_table TEXT NOT NULL,
    source_id UUID,
    error_type TEXT NOT NULL,
    error_details TEXT,
    original_data JSONB,
    suggested_fix JSONB,
    status TEXT DEFAULT 'pending',
    priority TEXT DEFAULT 'normal',
    reviewed_by UUID,
    reviewed_at TIMESTAMPTZ,
    review_notes TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_quarantined_status ON quarantined_data(status);
CREATE INDEX IF NOT EXISTS idx_quarantined_error ON quarantined_data(error_type);
CREATE INDEX IF NOT EXISTS idx_quarantined_source ON quarantined_data(source_table);


-- ============================================================================
-- 8. AI Fix Suggestions Table (AI 수정 제안)
-- ============================================================================

CREATE TABLE IF NOT EXISTS ai_fix_suggestions (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    quarantined_id UUID REFERENCES quarantined_data(id) ON DELETE CASCADE,
    suggestion_type TEXT NOT NULL,
    original_value TEXT,
    suggested_value TEXT,
    ai_model TEXT,
    ai_confidence NUMERIC(3,2),
    ai_reasoning TEXT,
    status TEXT DEFAULT 'pending',
    applied_by UUID,
    applied_at TIMESTAMPTZ,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

CREATE INDEX IF NOT EXISTS idx_ai_suggestions_status ON ai_fix_suggestions(status);
CREATE INDEX IF NOT EXISTS idx_ai_suggestions_quarantined ON ai_fix_suggestions(quarantined_id);


-- ============================================================================
-- 9. Live Session Stats (실시간 세션 통계 스냅샷)
-- ============================================================================

CREATE TABLE IF NOT EXISTS live_session_stats (
    id SERIAL PRIMARY KEY,
    snapshot_time TIMESTAMPTZ DEFAULT NOW(),
    total_active INT DEFAULT 0,
    by_stage JSONB,
    avg_duration_by_stage JSONB,
    bottleneck_stage TEXT,
    bottleneck_avg_seconds NUMERIC(10,2)
);

CREATE INDEX IF NOT EXISTS idx_live_stats_time ON live_session_stats(snapshot_time DESC);

-- 24시간 이전 데이터 자동 삭제 (옵션)
-- CREATE FUNCTION cleanup_old_stats() RETURNS void AS $$
-- DELETE FROM live_session_stats WHERE snapshot_time < NOW() - INTERVAL '24 hours';
-- $$ LANGUAGE SQL;


-- ============================================================================
-- 10. Data Health Metrics View (SKIP - 실제 테이블 스키마 확인 후 생성)
-- ============================================================================
-- 아래 VIEW는 실제 테이블의 컬럼명과 일치하지 않을 수 있으므로 주석 처리
-- 필요시 실제 스키마에 맞게 수정 후 실행

/*
CREATE OR REPLACE VIEW data_health_metrics AS
SELECT
    'antibody_library' AS source_table,
    COUNT(*) AS total_count,
    COUNT(*) FILTER (WHERE target_normalized IS NOT NULL) AS normalized_count,
    COUNT(*) FILTER (WHERE embedding IS NOT NULL) AS embedded_count,
    ROUND(COUNT(*) FILTER (WHERE target_normalized IS NOT NULL)::NUMERIC / NULLIF(COUNT(*), 0) * 100, 2) AS normalization_rate,
    ROUND(COUNT(*) FILTER (WHERE embedding IS NOT NULL)::NUMERIC / NULLIF(COUNT(*), 0) * 100, 2) AS embedding_rate
FROM antibody_library

UNION ALL

SELECT
    'commercial_reagents' AS source_table,
    COUNT(*) AS total_count,
    COUNT(*) FILTER (WHERE target_normalized IS NOT NULL) AS normalized_count,
    COUNT(*) FILTER (WHERE smiles IS NOT NULL AND smiles != '') AS with_smiles_count,
    ROUND(COUNT(*) FILTER (WHERE target_normalized IS NOT NULL)::NUMERIC / NULLIF(COUNT(*), 0) * 100, 2) AS normalization_rate,
    ROUND(COUNT(*) FILTER (WHERE smiles IS NOT NULL AND smiles != '')::NUMERIC / NULLIF(COUNT(*), 0) * 100, 2) AS smiles_rate
FROM commercial_reagents

UNION ALL

SELECT
    'golden_set_library' AS source_table,
    COUNT(*) AS total_count,
    COUNT(*) FILTER (WHERE target IS NOT NULL) AS normalized_count,
    COUNT(*) FILTER (WHERE orr_pct IS NOT NULL) AS with_clinical_count,
    ROUND(COUNT(*) FILTER (WHERE target IS NOT NULL)::NUMERIC / NULLIF(COUNT(*), 0) * 100, 2) AS normalization_rate,
    ROUND(COUNT(*) FILTER (WHERE orr_pct IS NOT NULL)::NUMERIC / NULLIF(COUNT(*), 0) * 100, 2) AS clinical_rate
FROM golden_set_library;
*/


-- ============================================================================
-- 11. Helper Functions
-- ============================================================================

-- 예산 사용량 업데이트 함수
CREATE OR REPLACE FUNCTION update_api_usage(
    p_provider TEXT,
    p_cost NUMERIC
)
RETURNS void
LANGUAGE plpgsql
AS $$
BEGIN
    UPDATE budget_configs
    SET current_usage_usd = current_usage_usd + p_cost,
        current_monthly_usd = current_monthly_usd + p_cost,
        updated_at = NOW()
    WHERE api_provider = p_provider;
END;
$$;


-- 일일 예산 리셋 함수
CREATE OR REPLACE FUNCTION reset_daily_budget()
RETURNS void
LANGUAGE plpgsql
AS $$
BEGIN
    UPDATE budget_configs
    SET current_usage_usd = 0,
        last_daily_reset_at = NOW()
    WHERE DATE(last_daily_reset_at) < CURRENT_DATE;
END;
$$;


-- 예산 초과 체크 함수
CREATE OR REPLACE FUNCTION check_budget_status(p_provider TEXT)
RETURNS TABLE (
    usage_pct NUMERIC,
    is_over_threshold BOOLEAN,
    is_over_limit BOOLEAN,
    should_fallback BOOLEAN
)
LANGUAGE plpgsql
AS $$
BEGIN
    RETURN QUERY
    SELECT
        ROUND((current_usage_usd / NULLIF(daily_limit_usd, 0)) * 100, 2) AS usage_pct,
        (current_usage_usd / NULLIF(daily_limit_usd, 0) * 100) >= alert_threshold_pct AS is_over_threshold,
        current_usage_usd >= daily_limit_usd AS is_over_limit,
        current_usage_usd >= daily_limit_usd AND auto_fallback_enabled AS should_fallback
    FROM budget_configs
    WHERE api_provider = p_provider;
END;
$$;


-- ============================================================================
-- 12. Permissions
-- ============================================================================

-- 권한은 RLS 정책에 따라 에러가 발생할 수 있음
-- 에러 발생 시 무시하고 계속 진행

DO $$
BEGIN
    GRANT SELECT, INSERT, UPDATE ON budget_configs TO authenticated;
EXCEPTION WHEN OTHERS THEN NULL;
END $$;

DO $$
BEGIN
    GRANT SELECT, INSERT ON api_usage_logs TO authenticated;
EXCEPTION WHEN OTHERS THEN NULL;
END $$;

DO $$
BEGIN
    GRANT SELECT, INSERT ON admin_audit_logs TO authenticated;
EXCEPTION WHEN OTHERS THEN NULL;
END $$;

DO $$
BEGIN
    GRANT SELECT, INSERT, UPDATE, DELETE ON target_synonyms TO authenticated;
EXCEPTION WHEN OTHERS THEN NULL;
END $$;

DO $$
BEGIN
    GRANT SELECT, INSERT, UPDATE ON quarantined_data TO authenticated;
EXCEPTION WHEN OTHERS THEN NULL;
END $$;

DO $$
BEGIN
    GRANT SELECT, INSERT, UPDATE ON ai_fix_suggestions TO authenticated;
EXCEPTION WHEN OTHERS THEN NULL;
END $$;

DO $$
BEGIN
    GRANT SELECT, INSERT ON live_session_stats TO authenticated;
EXCEPTION WHEN OTHERS THEN NULL;
END $$;

DO $$
BEGIN
    GRANT EXECUTE ON FUNCTION update_api_usage TO authenticated;
EXCEPTION WHEN OTHERS THEN NULL;
END $$;

DO $$
BEGIN
    GRANT EXECUTE ON FUNCTION reset_daily_budget TO authenticated;
EXCEPTION WHEN OTHERS THEN NULL;
END $$;

DO $$
BEGIN
    GRANT EXECUTE ON FUNCTION check_budget_status TO authenticated;
EXCEPTION WHEN OTHERS THEN NULL;
END $$;
