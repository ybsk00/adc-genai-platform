-- [Security Fix] 1. Enable RLS and add policies (Idempotent)
-- data_sync_logs
ALTER TABLE public.data_sync_logs ENABLE ROW LEVEL SECURITY;
DO $$ BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_policies WHERE tablename = 'data_sync_logs' AND policyname = 'Service Role Full Access') THEN
        CREATE POLICY "Service Role Full Access" ON public.data_sync_logs FOR ALL TO service_role USING (true) WITH CHECK (true);
    END IF;
    IF NOT EXISTS (SELECT 1 FROM pg_policies WHERE tablename = 'data_sync_logs' AND policyname = 'Authenticated Read Access') THEN
        CREATE POLICY "Authenticated Read Access" ON public.data_sync_logs FOR SELECT TO authenticated USING (true);
    END IF;
END $$;

-- transactions
ALTER TABLE public.transactions ENABLE ROW LEVEL SECURITY;
DO $$ BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_policies WHERE tablename = 'transactions' AND policyname = 'Service Role Full Access') THEN
        CREATE POLICY "Service Role Full Access" ON public.transactions FOR ALL TO service_role USING (true) WITH CHECK (true);
    END IF;
END $$;

-- kb_embeddings
ALTER TABLE public.kb_embeddings ENABLE ROW LEVEL SECURITY;
DO $$ BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_policies WHERE tablename = 'kb_embeddings' AND policyname = 'Service Role Full Access') THEN
        CREATE POLICY "Service Role Full Access" ON public.kb_embeddings FOR ALL TO service_role USING (true) WITH CHECK (true);
    END IF;
END $$;

-- commercial_reagents
ALTER TABLE public.commercial_reagents ENABLE ROW LEVEL SECURITY;
DO $$ BEGIN
    IF NOT EXISTS (SELECT 1 FROM pg_policies WHERE tablename = 'commercial_reagents' AND policyname = 'Service Role Full Access') THEN
        CREATE POLICY "Service Role Full Access" ON public.commercial_reagents FOR ALL TO service_role USING (true) WITH CHECK (true);
    END IF;
    IF NOT EXISTS (SELECT 1 FROM pg_policies WHERE tablename = 'commercial_reagents' AND policyname = 'Public Read Access') THEN
        CREATE POLICY "Public Read Access" ON public.commercial_reagents FOR SELECT TO anon, authenticated USING (true);
    END IF;
END $$;


-- [Security Fix] 2. Set Search Path for Functions (Correct Signatures)
DO $$
BEGIN
    -- approve_golden_set(uuid, uuid, text)
    IF EXISTS (SELECT 1 FROM pg_proc WHERE proname = 'approve_golden_set') THEN
        BEGIN
            EXECUTE 'ALTER FUNCTION public.approve_golden_set(uuid, uuid, text) SET search_path = public';
        EXCEPTION WHEN OTHERS THEN
            RAISE NOTICE 'Could not alter approve_golden_set: %', SQLERRM;
        END;
    END IF;

    -- handle_new_user()
    IF EXISTS (SELECT 1 FROM pg_proc WHERE proname = 'handle_new_user') THEN
        BEGIN
            EXECUTE 'ALTER FUNCTION public.handle_new_user() SET search_path = public';
        EXCEPTION WHEN OTHERS THEN
            RAISE NOTICE 'Could not alter handle_new_user: %', SQLERRM;
        END;
    END IF;

    -- update_updated_at()
    IF EXISTS (SELECT 1 FROM pg_proc WHERE proname = 'update_updated_at') THEN
        BEGIN
            EXECUTE 'ALTER FUNCTION public.update_updated_at() SET search_path = public';
        EXCEPTION WHEN OTHERS THEN
            RAISE NOTICE 'Could not alter update_updated_at: %', SQLERRM;
        END;
    END IF;

    -- match_golden_set(vector, float, int)
    IF EXISTS (SELECT 1 FROM pg_proc WHERE proname = 'match_golden_set') THEN
        BEGIN
            EXECUTE 'ALTER FUNCTION public.match_golden_set(vector, float, int) SET search_path = public';
        EXCEPTION WHEN OTHERS THEN
            RAISE NOTICE 'Could not alter match_golden_set: %', SQLERRM;
        END;
    END IF;
END
$$;
