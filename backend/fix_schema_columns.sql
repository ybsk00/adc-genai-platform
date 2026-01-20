-- Fix missing columns in golden_set_library
-- The table exists but might be missing 'raw_data' or other columns.

ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS raw_data JSONB DEFAULT '{}'::jsonb;

ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS properties JSONB DEFAULT '{}'::jsonb;

ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS enrichment_source TEXT;

ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS outcome_type TEXT;

ALTER TABLE public.golden_set_library 
ADD COLUMN IF NOT EXISTS failure_reason TEXT;

-- Re-apply permissions just in case
GRANT ALL ON public.golden_set_library TO service_role;
GRANT SELECT ON public.golden_set_library TO authenticated;
GRANT SELECT ON public.golden_set_library TO anon;
