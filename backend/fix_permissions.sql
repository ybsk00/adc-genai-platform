-- Fix permissions for golden_set_library
-- This script grants necessary table-level permissions and refreshes RLS policies.

-- 1. Grant Table Permissions
GRANT ALL ON public.golden_set_library TO service_role;
GRANT SELECT ON public.golden_set_library TO authenticated;
GRANT SELECT ON public.golden_set_library TO anon;

-- 2. Ensure RLS is enabled
ALTER TABLE public.golden_set_library ENABLE ROW LEVEL SECURITY;

-- 3. Refresh Policies
-- Service Role (Backend)
DROP POLICY IF EXISTS "Service Role Full Access" ON public.golden_set_library;
CREATE POLICY "Service Role Full Access" ON public.golden_set_library FOR ALL TO service_role USING (true) WITH CHECK (true);

-- Authenticated Users (Read)
DROP POLICY IF EXISTS "Authenticated Read Access" ON public.golden_set_library;
CREATE POLICY "Authenticated Read Access" ON public.golden_set_library FOR SELECT TO authenticated USING (true);

-- Admin Users (Write)
DROP POLICY IF EXISTS "Admin Full Access" ON public.golden_set_library;
CREATE POLICY "Admin Full Access" ON public.golden_set_library FOR ALL TO authenticated USING (true) WITH CHECK (true);
