-- [Fix] Create missing staging table for Data Sync
-- This table is used by scheduler.py to store drafted data from ClinicalTrials/PubMed

CREATE TABLE IF NOT EXISTS public.golden_set_library (
  id UUID DEFAULT gen_random_uuid() PRIMARY KEY,
  name TEXT NOT NULL,
  category TEXT,
  description TEXT,
  properties JSONB DEFAULT '{}'::jsonb,
  status TEXT DEFAULT 'draft', -- 'draft', 'approved', 'rejected'
  enrichment_source TEXT,
  raw_data JSONB DEFAULT '{}'::jsonb,
  outcome_type TEXT, -- 'Success', 'Failure', 'Ongoing', 'Unknown'
  failure_reason TEXT,
  created_at TIMESTAMPTZ DEFAULT NOW(),
  updated_at TIMESTAMPTZ DEFAULT NOW()
);

-- Enable RLS
ALTER TABLE public.golden_set_library ENABLE ROW LEVEL SECURITY;

-- Policies
DO $$ BEGIN
    -- Service Role (Backend) - Full Access
    IF NOT EXISTS (SELECT 1 FROM pg_policies WHERE tablename = 'golden_set_library' AND policyname = 'Service Role Full Access') THEN
        CREATE POLICY "Service Role Full Access" ON public.golden_set_library FOR ALL TO service_role USING (true) WITH CHECK (true);
    END IF;
    
    -- Authenticated Users (Frontend/Admin) - Read Access
    IF NOT EXISTS (SELECT 1 FROM pg_policies WHERE tablename = 'golden_set_library' AND policyname = 'Authenticated Read Access') THEN
        CREATE POLICY "Authenticated Read Access" ON public.golden_set_library FOR SELECT TO authenticated USING (true);
    END IF;
    
    -- Authenticated Users (Admin) - Write Access (for Approve/Reject)
    IF NOT EXISTS (SELECT 1 FROM pg_policies WHERE tablename = 'golden_set_library' AND policyname = 'Admin Full Access') THEN
        CREATE POLICY "Admin Full Access" ON public.golden_set_library FOR ALL TO authenticated USING (true) WITH CHECK (true);
    END IF;
END $$;

-- Indexes
CREATE INDEX IF NOT EXISTS idx_golden_set_library_status ON public.golden_set_library(status);
CREATE INDEX IF NOT EXISTS idx_golden_set_library_name ON public.golden_set_library(name);
