-- Phase 0: Closed-loop R&D Pipeline Schema

-- 1. design_runs table
-- Stores the configuration and status of a design run (Eng-Fit execution).
-- 'frozen_params' stores the snapshot of scoring parameters for auditability.
CREATE TABLE IF NOT EXISTS design_runs (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    project_id UUID REFERENCES projects(id),
    status TEXT CHECK (status IN ('draft', 'computing', 'completed', 'failed')) DEFAULT 'draft',
    frozen_params JSONB, -- Snapshot of scoring_params
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- 2. assay_results table
-- Stores the results of assays (virtual or wet-lab) linked to a design run.
-- 'acceptance_criteria' stores the rule used to judge success.
CREATE TABLE IF NOT EXISTS assay_results (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    run_id UUID REFERENCES design_runs(id),
    molecule_id UUID REFERENCES golden_set_library(id), -- Linking to ADC candidates
    assay_type TEXT,
    raw_data JSONB,
    is_success BOOLEAN,
    acceptance_criteria JSONB,
    confidence_score FLOAT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- 3. component_catalog table
-- Stores individual components (Antibody, Payload, Linker) and their properties.
-- Added locking columns for Async-Precompute.
CREATE TABLE IF NOT EXISTS component_catalog (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT,
    smiles TEXT,
    category TEXT, -- 'antibody', 'payload', 'linker'
    properties JSONB,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Add locking columns if they don't exist (Idempotent)
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'component_catalog' AND column_name = 'lock_status') THEN
        ALTER TABLE component_catalog ADD COLUMN lock_status TEXT DEFAULT 'available'; -- 'available', 'computing'
    END IF;

    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'component_catalog' AND column_name = 'lock_holder') THEN
        ALTER TABLE component_catalog ADD COLUMN lock_holder UUID; -- Worker ID
    END IF;

    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'component_catalog' AND column_name = 'estimated_wait_time') THEN
        ALTER TABLE component_catalog ADD COLUMN estimated_wait_time INTEGER; -- Seconds
    END IF;
END $$;
