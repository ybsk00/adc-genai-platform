-- [Migration] Final Advanced Schema Update for golden_set_library
-- This migration ensures all necessary columns for efficacy extraction and ADC metrics are present.

-- 1. Clinical Figures (Efficacy & Safety)
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS orr_pct NUMERIC,
ADD COLUMN IF NOT EXISTS os_months NUMERIC,
ADD COLUMN IF NOT EXISTS pfs_months NUMERIC,
ADD COLUMN IF NOT EXISTS dor_months NUMERIC,
ADD COLUMN IF NOT EXISTS patient_count INTEGER,
ADD COLUMN IF NOT EXISTS adverse_events_grade3_pct NUMERIC;

-- 2. ADC Specifics & Bio Metrics
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS dar NUMERIC,
ADD COLUMN IF NOT EXISTS binding_affinity TEXT,
ADD COLUMN IF NOT EXISTS isotype TEXT,
ADD COLUMN IF NOT EXISTS host_species TEXT,
ADD COLUMN IF NOT EXISTS antibody_format TEXT,
ADD COLUMN IF NOT EXISTS molecular_weight TEXT;

-- 3. Chemical Structure
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS payload_smiles TEXT,
ADD COLUMN IF NOT EXISTS linker_smiles TEXT,
ADD COLUMN IF NOT EXISTS full_smiles TEXT,
ADD COLUMN IF NOT EXISTS linker_type TEXT;

-- 4. Bio Targets
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS target_1 TEXT,
ADD COLUMN IF NOT EXISTS target_2 TEXT,
ADD COLUMN IF NOT EXISTS gene_id TEXT,
ADD COLUMN IF NOT EXISTS uniprot_id TEXT;

-- 5. Quality Control & Metadata
ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS confidence_score FLOAT DEFAULT 0.0,
ADD COLUMN IF NOT EXISTS review_required BOOLEAN DEFAULT FALSE,
ADD COLUMN IF NOT EXISTS ai_refined BOOLEAN DEFAULT FALSE,
ADD COLUMN IF NOT EXISTS rag_status TEXT DEFAULT 'pending',
ADD COLUMN IF NOT EXISTS processing_error TEXT;

-- 6. Indexes for Optimization
CREATE INDEX IF NOT EXISTS idx_gsl_target_1 ON public.golden_set_library(target_1);
CREATE INDEX IF NOT EXISTS idx_gsl_ai_refined ON public.golden_set_library(ai_refined);
CREATE INDEX IF NOT EXISTS idx_gsl_review_required ON public.golden_set_library(review_required);
CREATE INDEX IF NOT EXISTS idx_gsl_outcome_type ON public.golden_set_library(outcome_type);
