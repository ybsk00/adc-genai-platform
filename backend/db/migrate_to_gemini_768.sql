-- Migration: 1536 -> 768 Dimensions (Gemini text-embedding-004)

-- 1. Drop existing functions that depend on vector(1536)
DROP FUNCTION IF EXISTS match_golden_set(vector, float, int);
DROP FUNCTION IF EXISTS match_golden_set(vector, double precision, integer);

-- 2. Update golden_set_embeddings table
ALTER TABLE golden_set_embeddings DROP COLUMN IF EXISTS embedding;
ALTER TABLE golden_set_embeddings ADD COLUMN embedding vector(768);

-- 3. Update kb_embeddings table (if exists)
DO $$
BEGIN
    IF EXISTS (SELECT 1 FROM information_schema.tables WHERE table_name = 'kb_embeddings') THEN
        ALTER TABLE kb_embeddings DROP COLUMN IF EXISTS embedding;
        ALTER TABLE kb_embeddings ADD COLUMN embedding vector(768);
    END IF;
END $$;

-- 4. Update commercial_reagents table
DO $$
BEGIN
    IF EXISTS (SELECT 1 FROM information_schema.tables WHERE table_name = 'commercial_reagents') THEN
        ALTER TABLE commercial_reagents DROP COLUMN IF EXISTS embedding;
        ALTER TABLE commercial_reagents ADD COLUMN embedding vector(768);
        
        -- Create index for commercial_reagents
        CREATE INDEX IF NOT EXISTS idx_commercial_reagents_embedding ON commercial_reagents USING hnsw (embedding vector_cosine_ops);
    END IF;
END $$;

-- 5. Recreate match_golden_set function with vector(768)
CREATE OR REPLACE FUNCTION match_golden_set (
  query_embedding vector(768),
  match_threshold float,
  match_count int
)
RETURNS TABLE (
  id uuid,
  content text,
  similarity float,
  metadata jsonb
)
LANGUAGE plpgsql
AS $$
BEGIN
  RETURN QUERY
  SELECT
    g.id,
    e.chunk_content AS content,
    1 - (e.embedding <=> query_embedding) AS similarity,
    g.properties AS metadata
  FROM golden_set_embeddings e
  JOIN golden_set_library g ON e.source_id = g.id
  WHERE 1 - (e.embedding <=> query_embedding) > match_threshold
  AND g.status = 'approved'
  ORDER BY similarity DESC
  LIMIT match_count;
END;
$$;
