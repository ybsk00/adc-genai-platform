-- Fix embedding dimension mismatch (1536 -> 768)

-- 1. Drop existing column (if exists)
ALTER TABLE commercial_reagents DROP COLUMN IF EXISTS embedding;

-- 2. Add new column with 768 dimensions
-- Note: Requires pgvector extension (usually enabled in Supabase)
ALTER TABLE commercial_reagents ADD COLUMN embedding vector(768);

-- 3. Create index for faster similarity search
CREATE INDEX IF NOT EXISTS idx_commercial_reagents_embedding ON commercial_reagents USING hnsw (embedding vector_cosine_ops);
