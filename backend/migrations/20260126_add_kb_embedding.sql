
-- knowledge_base 테이블에 Gemini 768차원 임베딩 컬럼 추가
ALTER TABLE public.knowledge_base ADD COLUMN IF NOT EXISTS embedding vector(768);

-- 벡터 검색 성능 향상을 위한 인덱스 추가 (HNSW)
CREATE INDEX IF NOT EXISTS idx_knowledge_base_embedding 
ON public.knowledge_base USING hnsw (embedding vector_cosine_ops);

-- RLS 정책 확인 (기존 정책 유지)
COMMENT ON COLUMN public.knowledge_base.embedding IS 'Gemini 2.0 Flash Generated 768-dim Embedding';
