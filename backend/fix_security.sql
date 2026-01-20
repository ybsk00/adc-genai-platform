-- [Security Fix] 1. Enable RLS and add policies
-- data_sync_logs
ALTER TABLE public.data_sync_logs ENABLE ROW LEVEL SECURITY;
CREATE POLICY "Service Role Full Access" ON public.data_sync_logs FOR ALL TO service_role USING (true) WITH CHECK (true);
CREATE POLICY "Authenticated Read Access" ON public.data_sync_logs FOR SELECT TO authenticated USING (true);

-- transactions
ALTER TABLE public.transactions ENABLE ROW LEVEL SECURITY;
CREATE POLICY "Service Role Full Access" ON public.transactions FOR ALL TO service_role USING (true) WITH CHECK (true);
-- Note: Add user-specific policies later if needed

-- kb_embeddings
ALTER TABLE public.kb_embeddings ENABLE ROW LEVEL SECURITY;
CREATE POLICY "Service Role Full Access" ON public.kb_embeddings FOR ALL TO service_role USING (true) WITH CHECK (true);

-- commercial_reagents
ALTER TABLE public.commercial_reagents ENABLE ROW LEVEL SECURITY;
CREATE POLICY "Service Role Full Access" ON public.commercial_reagents FOR ALL TO service_role USING (true) WITH CHECK (true);
CREATE POLICY "Public Read Access" ON public.commercial_reagents FOR SELECT TO anon, authenticated USING (true);


-- [Security Fix] 2. Set Search Path for Functions (to prevent search_path hijacking)
ALTER FUNCTION public.approve_golden_set(uuid) SET search_path = public;
ALTER FUNCTION public.handle_new_user() SET search_path = public;
ALTER FUNCTION public.update_updated_at() SET search_path = public;
-- Note: match_golden_set signature might vary, assuming generic or checking existence
DO $$
BEGIN
    IF EXISTS (SELECT 1 FROM pg_proc WHERE proname = 'match_golden_set') THEN
        EXECUTE 'ALTER FUNCTION public.match_golden_set SET search_path = public';
    END IF;
END
$$;
