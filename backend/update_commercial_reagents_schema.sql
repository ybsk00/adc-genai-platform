-- Add new columns to commercial_reagents table for AI refinement and RAG
-- Safe to run multiple times (IF NOT EXISTS)

DO $$
BEGIN
    -- ai_refined
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'commercial_reagents' AND column_name = 'ai_refined') THEN
        ALTER TABLE commercial_reagents ADD COLUMN ai_refined BOOLEAN DEFAULT FALSE;
    END IF;

    -- target
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'commercial_reagents' AND column_name = 'target') THEN
        ALTER TABLE commercial_reagents ADD COLUMN target TEXT;
    END IF;

    -- properties (JSONB) - might already exist, but ensuring
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'commercial_reagents' AND column_name = 'properties') THEN
        ALTER TABLE commercial_reagents ADD COLUMN properties JSONB DEFAULT '{}'::jsonb;
    END IF;

    -- summary
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'commercial_reagents' AND column_name = 'summary') THEN
        ALTER TABLE commercial_reagents ADD COLUMN summary TEXT;
    END IF;

    -- source_name
    IF NOT EXISTS (SELECT 1 FROM information_schema.columns WHERE table_name = 'commercial_reagents' AND column_name = 'source_name') THEN
        ALTER TABLE commercial_reagents ADD COLUMN source_name TEXT;
    END IF;

    -- Create index for ai_refined for faster querying
    CREATE INDEX IF NOT EXISTS idx_commercial_reagents_ai_refined ON commercial_reagents(ai_refined);
    
    -- Create index for source_name
    CREATE INDEX IF NOT EXISTS idx_commercial_reagents_source_name ON commercial_reagents(source_name);

END $$;
