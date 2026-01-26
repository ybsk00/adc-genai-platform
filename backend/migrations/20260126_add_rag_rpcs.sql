-- 1. Knowledge Base Search Function
CREATE OR REPLACE FUNCTION match_knowledge_base (
  query_embedding vector(768),
  match_threshold float,
  match_count int
)
RETURNS TABLE (
  id bigint,
  content text,
  similarity float,
  metadata jsonb
)
LANGUAGE plpgsql
AS $$
BEGIN
  RETURN QUERY
  SELECT
    kb.id,
    kb.content,
    1 - (kb.embedding <=> query_embedding) as similarity,
    jsonb_build_object(
      'title', kb.title,
      'source_type', kb.source_type,
      'summary', kb.summary
    ) as metadata
  FROM knowledge_base kb
  WHERE 1 - (kb.embedding <=> query_embedding) > match_threshold
  ORDER BY similarity DESC
  LIMIT match_count;
END;
$$;

-- 2. Antibody Library Search Function
CREATE OR REPLACE FUNCTION match_antibody_library (
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
    al.id,
    al.summary as content, -- Use summary as the main content for RAG
    1 - (al.embedding <=> query_embedding) as similarity,
    jsonb_build_object(
      'product_name', al.product_name,
      'cat_no', al.cat_no,
      'host_species', al.host_species,
      'isotype', al.isotype
    ) as metadata
  FROM antibody_library al
  WHERE 1 - (al.embedding <=> query_embedding) > match_threshold
  ORDER BY similarity DESC
  LIMIT match_count;
END;
$$;
