-- [Migration] Add Advanced Columns to golden_set_library
-- Clinical Figures
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS dor_months NUMERIC,
ADD COLUMN IF NOT EXISTS patient_count INTEGER,
ADD COLUMN IF NOT EXISTS adverse_events_grade3_pct NUMERIC;

-- ADC Specifics
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS dar NUMERIC;

-- Chemical Structure
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS payload_smiles TEXT,
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS linker_smiles TEXT,
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS full_smiles TEXT,
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS linker_type TEXT;

-- Bio Targets (Bispecifics & Details)
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS target_1 TEXT,
ADD COLUMN IF NOT EXISTS target_2 TEXT,
ADD COLUMN IF NOT EXISTS gene_id TEXT,
ADD COLUMN IF NOT EXISTS uniprot_id TEXT,
ADD COLUMN IF NOT EXISTS antibody_format TEXT,
ADD COLUMN IF NOT EXISTS molecular_weight TEXT;

-- Quality Control
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS confidence_score FLOAT DEFAULT 0.0,
ADD COLUMN IF NOT EXISTS review_required BOOLEAN DEFAULT FALSE;

-- Indexes for frequent queries
CREATE INDEX IF NOT EXISTS idx_gsl_target_1 ON public.golden_set_library(target_1);
CREATE INDEX IF NOT EXISTS idx_gsl_review_required ON public.golden_set_library(review_required);
