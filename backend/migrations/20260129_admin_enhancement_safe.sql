-- ============================================================================
-- Migration: Admin Dashboard Enhancement v2.2+ (Safe Version)
-- Date: 2026-01-29
-- 기존 system_config 테이블이 있는 경우를 위한 안전한 마이그레이션
--
-- 실행 순서: 섹션별로 하나씩 실행하세요
-- ============================================================================


-- ============================================================================
-- STEP 1: 기존 system_config에 설정값 추가
-- ============================================================================
INSERT INTO system_config (key, value) VALUES
    ('USE_NIM_API', '"true"')
ON CONFLICT (key) DO NOTHING;

INSERT INTO system_config (key, value) VALUES
    ('FALLBACK_ENABLED', '"true"')
ON CONFLICT (key) DO NOTHING;

INSERT INTO system_config (key, value) VALUES
    ('EMERGENCY_STOP', '"false"')
ON CONFLICT (key) DO NOTHING;


-- ============================================================================
-- STEP 2: Budget Configuration Table (API 비용 관리)
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

INSERT INTO budget_configs (api_provider, daily_limit_usd, monthly_limit_usd) VALUES
    ('nvidia_nim', 500.00, 10000.00),
    ('gemini', 200.00, 4000.00),
    ('alphafold', 300.00, 6000.00)
ON CONFLICT (api_provider) DO NOTHING;


-- ============================================================================
-- STEP 3: API Usage Logs Table (API 호출 로그)
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
-- STEP 4: Admin Audit Logs Table (관리자 감사 로그)
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
-- STEP 5: Target Synonyms Table (타겟 동의어 사전)
-- ============================================================================
CREATE TABLE IF NOT EXISTS target_synonyms (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    canonical_name TEXT NOT NULL UNIQUE,
    synonyms TEXT[] NOT NULL,
    category TEXT,
    description TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

INSERT INTO target_synonyms (canonical_name, synonyms, category) VALUES
    ('HER2', ARRAY['ERBB2', 'neu', 'CD340', 'HER-2', 'HER2/neu'], 'receptor'),
    ('EGFR', ARRAY['HER1', 'ERBB1', 'ErbB-1'], 'receptor'),
    ('PD-1', ARRAY['PD1', 'PDCD1', 'CD279'], 'checkpoint'),
    ('PD-L1', ARRAY['PDL1', 'CD274', 'B7-H1'], 'checkpoint'),
    ('TROP2', ARRAY['TROP-2', 'TACSTD2', 'EGP-1'], 'antigen'),
    ('CD19', ARRAY['B4', 'CVID3'], 'antigen'),
    ('CD33', ARRAY['SIGLEC3', 'gp67'], 'antigen'),
    ('BCMA', ARRAY['TNFRSF17', 'CD269'], 'antigen'),
    ('VEGF', ARRAY['VEGFA', 'VPF'], 'growth_factor'),
    ('CD3', ARRAY['T3'], 'receptor')
ON CONFLICT (canonical_name) DO NOTHING;


-- ============================================================================
-- STEP 6: Quarantined Data Table (격리 데이터)
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
-- STEP 7: AI Fix Suggestions Table (AI 수정 제안)
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
-- STEP 8: Live Session Stats (실시간 세션 통계 스냅샷)
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


-- ============================================================================
-- STEP 9: Permissions (RLS가 활성화된 경우)
-- ============================================================================
-- authenticated 역할에 권한 부여
-- 에러가 나면 무시해도 됨 (RLS가 비활성화되어 있을 수 있음)

GRANT SELECT, INSERT, UPDATE ON budget_configs TO authenticated;
GRANT SELECT, INSERT ON api_usage_logs TO authenticated;
GRANT SELECT, INSERT ON admin_audit_logs TO authenticated;
GRANT SELECT, INSERT, UPDATE, DELETE ON target_synonyms TO authenticated;
GRANT SELECT, INSERT, UPDATE ON quarantined_data TO authenticated;
GRANT SELECT, INSERT, UPDATE ON ai_fix_suggestions TO authenticated;
GRANT SELECT, INSERT ON live_session_stats TO authenticated;


-- ============================================================================
-- 완료!
-- ============================================================================
-- 모든 테이블이 생성되었습니다.
-- Admin 대시보드에서 Engine Control, Refinement Hub, Audit 메뉴를 사용할 수 있습니다.
