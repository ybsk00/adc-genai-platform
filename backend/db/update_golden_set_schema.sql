-- Add Clinical Metrics
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS orr_pct numeric;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS os_months numeric;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS pfs_months numeric;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS dor_months numeric;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS patient_count integer;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS adverse_events_grade3_pct numeric;

-- Add ADC Specific Metrics
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS dar numeric;

-- Add Chemical Structure Fields
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS payload_smiles text;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS linker_smiles text;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS full_smiles text;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS linker_type text;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS antibody_format text;

-- Add Biological Targets (Bispecific Support)
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS target_1 text;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS target_2 text;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS gene_id text;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS uniprot_id text;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS binding_affinity text;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS isotype text;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS molecular_weight numeric;

-- Add Quality Control Fields
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS confidence_score numeric DEFAULT 0.0;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS review_required boolean DEFAULT false;

-- Create Indexes for frequent search columns
CREATE INDEX IF NOT EXISTS idx_golden_set_target_1 ON golden_set_library (target_1);
CREATE INDEX IF NOT EXISTS idx_golden_set_target_2 ON golden_set_library (target_2);
CREATE INDEX IF NOT EXISTS idx_golden_set_review_required ON golden_set_library (review_required);
