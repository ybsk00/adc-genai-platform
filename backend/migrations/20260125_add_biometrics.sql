-- [Migration] Add Bio Metrics columns to golden_set_library table
-- Run this in your Supabase SQL Editor

ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS binding_affinity TEXT,
ADD COLUMN IF NOT EXISTS isotype TEXT,
ADD COLUMN IF NOT EXISTS host_species TEXT,
ADD COLUMN IF NOT EXISTS orr_pct TEXT,
ADD COLUMN IF NOT EXISTS os_months TEXT,
ADD COLUMN IF NOT EXISTS pfs_months TEXT,
ADD COLUMN IF NOT EXISTS is_manual_override BOOLEAN DEFAULT FALSE;

-- Optional: Re-index for search optimization
CREATE INDEX IF NOT EXISTS idx_gsl_is_manual_override ON public.golden_set_library(is_manual_override);
