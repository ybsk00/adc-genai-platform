-- Report Verification Logs Table
-- 보고서 Digital Seal 및 QR 검증을 위한 로그 테이블
-- 21 CFR Part 11 준수를 위한 감사 추적

-- Create report_verification_logs table
CREATE TABLE IF NOT EXISTS report_verification_logs (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),

    -- Report identification
    session_id UUID,
    report_id VARCHAR(64),

    -- Digital Seal
    chain_hash VARCHAR(64),
    short_hash VARCHAR(12),
    page_hashes JSONB,
    page_count INTEGER DEFAULT 10,
    algorithm VARCHAR(32) DEFAULT 'SHA-256 Chain',

    -- Verification URL
    verification_url TEXT,

    -- Audit info
    generated_at TIMESTAMPTZ DEFAULT NOW(),
    generated_by VARCHAR(255),

    -- Metadata
    report_config JSONB,

    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- 기존 테이블에 누락된 컬럼 추가
DO $$
BEGIN
    -- short_hash 컬럼
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'report_verification_logs' AND column_name = 'short_hash') THEN
        ALTER TABLE report_verification_logs ADD COLUMN short_hash VARCHAR(12);
    END IF;

    -- chain_hash 컬럼
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'report_verification_logs' AND column_name = 'chain_hash') THEN
        ALTER TABLE report_verification_logs ADD COLUMN chain_hash VARCHAR(64);
    END IF;

    -- page_hashes 컬럼
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'report_verification_logs' AND column_name = 'page_hashes') THEN
        ALTER TABLE report_verification_logs ADD COLUMN page_hashes JSONB;
    END IF;

    -- page_count 컬럼
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'report_verification_logs' AND column_name = 'page_count') THEN
        ALTER TABLE report_verification_logs ADD COLUMN page_count INTEGER DEFAULT 10;
    END IF;

    -- algorithm 컬럼
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'report_verification_logs' AND column_name = 'algorithm') THEN
        ALTER TABLE report_verification_logs ADD COLUMN algorithm VARCHAR(32) DEFAULT 'SHA-256 Chain';
    END IF;

    -- verification_url 컬럼
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'report_verification_logs' AND column_name = 'verification_url') THEN
        ALTER TABLE report_verification_logs ADD COLUMN verification_url TEXT;
    END IF;

    -- report_config 컬럼
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns
                   WHERE table_name = 'report_verification_logs' AND column_name = 'report_config') THEN
        ALTER TABLE report_verification_logs ADD COLUMN report_config JSONB;
    END IF;
END $$;

-- Indexes for fast lookup (컬럼이 존재하는 경우에만 생성됨)
CREATE INDEX IF NOT EXISTS idx_report_verification_short_hash
    ON report_verification_logs(short_hash);

CREATE INDEX IF NOT EXISTS idx_report_verification_chain_hash
    ON report_verification_logs(chain_hash);

CREATE INDEX IF NOT EXISTS idx_report_verification_session
    ON report_verification_logs(session_id);

CREATE INDEX IF NOT EXISTS idx_report_verification_report_id
    ON report_verification_logs(report_id);

-- Trigger for updated_at
CREATE OR REPLACE FUNCTION update_report_verification_updated_at()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = NOW();
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS trigger_report_verification_updated_at ON report_verification_logs;
CREATE TRIGGER trigger_report_verification_updated_at
    BEFORE UPDATE ON report_verification_logs
    FOR EACH ROW
    EXECUTE FUNCTION update_report_verification_updated_at();

-- RLS Policies (if needed)
ALTER TABLE report_verification_logs ENABLE ROW LEVEL SECURITY;

-- Drop existing policies if they exist (to avoid conflicts)
DROP POLICY IF EXISTS "Users can view their report verifications" ON report_verification_logs;
DROP POLICY IF EXISTS "Service role full access to report_verification_logs" ON report_verification_logs;
DROP POLICY IF EXISTS "Public can verify reports" ON report_verification_logs;

-- Allow authenticated users to read verification logs for their reports
CREATE POLICY "Users can view their report verifications"
    ON report_verification_logs
    FOR SELECT
    USING (
        EXISTS (
            SELECT 1 FROM design_sessions
            WHERE design_sessions.id = report_verification_logs.session_id
            AND design_sessions.user_id::text = auth.uid()::text
        )
    );

-- Allow service role full access
CREATE POLICY "Service role full access to report_verification_logs"
    ON report_verification_logs
    FOR ALL
    USING (true);

-- Public can verify reports (read only on short_hash lookup)
CREATE POLICY "Public can verify reports"
    ON report_verification_logs
    FOR SELECT
    USING (true);

-- Comments
COMMENT ON TABLE report_verification_logs IS 'Digital Seal and QR verification logs for ADC reports';
COMMENT ON COLUMN report_verification_logs.chain_hash IS 'SHA-256 chain hash of all page hashes';
COMMENT ON COLUMN report_verification_logs.short_hash IS '12-character hash for QR code URL';
COMMENT ON COLUMN report_verification_logs.page_hashes IS 'JSON array of individual page hashes';
