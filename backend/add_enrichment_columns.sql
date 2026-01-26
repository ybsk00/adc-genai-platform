-- Add missing columns for Clinical Data Enrichment
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS patient_count INTEGER;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS dar FLOAT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS adverse_events_grade3_pct FLOAT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS target_symbol TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS target_1 TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS target_2 TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS gene_id TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS uniprot_id TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS antibody_format TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS binding_affinity TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS isotype TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS host_species TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS molecular_weight TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS dor_months FLOAT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS review_required BOOLEAN DEFAULT FALSE;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS payload_smiles TEXT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS confidence_score FLOAT;
ALTER TABLE golden_set_library ADD COLUMN IF NOT EXISTS relevance_score FLOAT;
