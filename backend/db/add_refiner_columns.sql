-- AI Refiner 고도화를 위한 컬럼 추가
-- 2026-01-21

DO $$
BEGIN
    -- 1. relevance_score (ADC 관련성 점수)
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'golden_set_library' AND column_name = 'relevance_score') THEN
        ALTER TABLE public.golden_set_library ADD COLUMN relevance_score FLOAT;
    END IF;

    -- 2. rag_status (처리 상태)
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'golden_set_library' AND column_name = 'rag_status') THEN
        ALTER TABLE public.golden_set_library ADD COLUMN rag_status TEXT DEFAULT 'pending';
    END IF;

    -- 3. molecular_weight (분자량)
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'golden_set_library' AND column_name = 'molecular_weight') THEN
        ALTER TABLE public.golden_set_library ADD COLUMN molecular_weight FLOAT;
    END IF;

    -- 4. canonical_smiles (표준 SMILES)
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'golden_set_library' AND column_name = 'canonical_smiles') THEN
        ALTER TABLE public.golden_set_library ADD COLUMN canonical_smiles TEXT;
    END IF;

    -- 5. ai_refined (기존 컬럼이지만 확실히 하기 위해)
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'golden_set_library' AND column_name = 'ai_refined') THEN
        ALTER TABLE public.golden_set_library ADD COLUMN ai_refined BOOLEAN DEFAULT FALSE;
    END IF;

    -- 6. processing_error (기존 컬럼이지만 확실히 하기 위해)
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'golden_set_library' AND column_name = 'processing_error') THEN
        ALTER TABLE public.golden_set_library ADD COLUMN processing_error TEXT;
    END IF;
END $$;
