-- [1] Target Master (타겟 정보)
CREATE TABLE IF NOT EXISTS public.target_master (
    id UUID NOT NULL DEFAULT gen_random_uuid(),
    target_symbol TEXT NOT NULL UNIQUE, -- e.g., 'HER2', 'CD3'
    gene_id TEXT,
    uniprot_id TEXT,
    alt_names TEXT[], -- Alternative Names (Array)
    target_type TEXT, -- e.g., 'Tumor Antigen', 'Immune Cell'
    created_at TIMESTAMP WITH TIME ZONE DEFAULT now(),
    CONSTRAINT target_master_pkey PRIMARY KEY (id)
);

-- [2] Antibody Library (항체 라이브러리 - 상세 스펙)
CREATE TABLE IF NOT EXISTS public.antibody_library (
    id UUID NOT NULL DEFAULT gen_random_uuid(),
    product_name TEXT NOT NULL,
    cat_no TEXT NOT NULL UNIQUE, -- Creative Biolabs Cat No.
    target_id UUID REFERENCES public.target_master(id),
    host_species TEXT,
    isotype TEXT,
    clone_id TEXT,
    conjugate TEXT,
    application TEXT[], -- e.g., ['Flow Cyt', 'WB']
    full_spec JSONB, -- 기타 상세 스펙 (Reactivity, Formulation 등)
    source_url TEXT,
    source_name TEXT DEFAULT 'Creative Biolabs',
    
    -- AI/Vector Search
    embedding vector(768), -- Gemini 768-dim
    summary TEXT,
    
    crawled_at TIMESTAMP WITH TIME ZONE DEFAULT now(),
    CONSTRAINT antibody_library_pkey PRIMARY KEY (id)
);

-- [3] Application Map (이중항체 등 복합 타겟 매핑용)
CREATE TABLE IF NOT EXISTS public.application_map (
    id UUID NOT NULL DEFAULT gen_random_uuid(),
    antibody_id UUID REFERENCES public.antibody_library(id) ON DELETE CASCADE,
    target_symbol TEXT,
    interaction_type TEXT, -- e.g., 'Primary Binder', 'Bispecific Arm 1'
    CONSTRAINT application_map_pkey PRIMARY KEY (id)
);

-- Indexes for Speed
CREATE INDEX IF NOT EXISTS idx_antibody_library_embedding ON public.antibody_library USING hnsw (embedding vector_cosine_ops);
CREATE INDEX IF NOT EXISTS idx_antibody_library_cat_no ON public.antibody_library(cat_no);
CREATE INDEX IF NOT EXISTS idx_target_master_symbol ON public.target_master(target_symbol);
