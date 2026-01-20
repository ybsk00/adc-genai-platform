-- 데이터 동기화 로그 테이블 생성
CREATE TABLE IF NOT EXISTS public.data_sync_logs (
    id BIGINT PRIMARY KEY GENERATED ALWAYS AS IDENTITY,
    source_id TEXT NOT NULL, -- 'clinical_trials', 'pubmed', 'goldenset' 등
    status TEXT NOT NULL,    -- 'success', 'failed', 'partial'
    records_synced INT DEFAULT 0,
    records_drafted INT DEFAULT 0,
    error_message TEXT,
    created_at TIMESTAMPTZ DEFAULT NOW()
);

-- RLS 설정 (관리자만 조회 가능하도록 설정할 수 있으나, 일단은 단순화)
ALTER TABLE public.data_sync_logs ENABLE ROW LEVEL SECURITY;
CREATE POLICY "Allow all for now" ON public.data_sync_logs FOR ALL USING (true);
