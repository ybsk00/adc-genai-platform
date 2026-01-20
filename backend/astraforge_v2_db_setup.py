print("Please run the following SQL in Supabase SQL Editor to upgrade to AstraForge 2.0:")
print("""
-- 1. 상용 시약 테이블 생성
CREATE TABLE IF NOT EXISTS commercial_reagents (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL,
    category TEXT, -- antibody, linker, payload, linker_payload
    description TEXT,
    supplier TEXT,
    product_url TEXT,
    catalog_no TEXT UNIQUE,
    price_data JSONB,
    properties JSONB DEFAULT '{}'::jsonb, -- uniprot_id, target 등 포함
    created_at TIMESTAMPTZ DEFAULT NOW(),
    updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- 2. 골든셋 라이브러리 필드 확장
ALTER TABLE golden_set_library 
ADD COLUMN IF NOT EXISTS outcome_type TEXT DEFAULT 'Unknown', -- Success, Failure, Ongoing
ADD COLUMN IF NOT EXISTS failure_reason TEXT,
ADD COLUMN IF NOT EXISTS is_ai_extracted BOOLEAN DEFAULT FALSE;

-- 3. RLS 정책 (필요시)
ALTER TABLE commercial_reagents ENABLE ROW LEVEL SECURITY;
CREATE POLICY "Allow all for authenticated users" ON commercial_reagents 
    FOR ALL TO authenticated USING (true) WITH CHECK (true);
""")
