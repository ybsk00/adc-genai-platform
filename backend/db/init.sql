-- 1. Enable Vector Extension (RAG 필수)
create extension if not exists vector;

-- 2. Profiles Table (사용자 정보)
create table if not exists public.profiles (
  id uuid references auth.users on delete cascade not null primary key,
  email text unique not null,
  full_name text,
  org_name text,
  credits int default 10, -- 가입 시 무료 크레딧 10 제공
  created_at timestamptz default now(),
  updated_at timestamptz default now()
);

-- 3. Projects Table (시뮬레이션 작업)
create table if not exists public.projects (
  id uuid default gen_random_uuid() primary key,
  user_id uuid references public.profiles(id) not null,
  title text not null,
  target_protein text, -- e.g. "LIV-1"
  payload_type text,   -- e.g. "MMAE"
  status text default 'queued', -- queued, processing, completed, failed
  result_summary jsonb, -- AI 결과 요약
  report_url text,      -- PDF 다운로드 링크
  created_at timestamptz default now()
);

-- 4. Golden Set Library (약물/임상 데이터)
-- Human-in-the-Loop 적용: status 컬럼 포함
create table if not exists public.golden_set_library (
  id uuid default gen_random_uuid() primary key,
  name text not null, -- 약물명 or 타겟명
  category text,      -- antibody, payload, linker, clinical_trial
  description text,   -- RAG 검색용 텍스트
  properties jsonb,   -- 상세 스펙 (IC50, Phase 등)
  
  -- 검증 프로세스용 컬럼
  status text default 'draft' check (status in ('draft', 'approved', 'rejected')),
  enrichment_source text, -- 'perplexity', 'clinical_trials'
  reviewer_note text,     -- 관리자 반려 사유
  
  created_at timestamptz default now()
);

-- 5. Golden Set Embeddings (벡터 저장소)
create table if not exists public.golden_set_embeddings (
  id uuid default gen_random_uuid() primary key,
  source_id uuid references public.golden_set_library(id) on delete cascade,
  chunk_content text,       -- 실제 임베딩된 텍스트 조각
  embedding vector(1536)    -- OpenAI Embedding Dimension
);

-- 6. Transactions (결제 및 크레딧 로그)
create table if not exists public.transactions (
  id uuid default gen_random_uuid() primary key,
  user_id uuid references public.profiles(id) not null,
  amount int not null,      -- 변동량 (+500, -10)
  type text not null,       -- purchase, usage, bonus
  description text,
  created_at timestamptz default now()
);

-- 7. RAG Search Function (벡터 검색 함수)
-- 기존 함수가 있으면 리턴 타입 변경이 불가능하므로 삭제 후 재생성
DROP FUNCTION IF EXISTS match_golden_set(vector, float, int);
DROP FUNCTION IF EXISTS match_golden_set(vector, double precision, integer);

create or replace function match_golden_set (
  query_embedding vector(1536),
  match_threshold float,
  match_count int
)
returns table (
  id uuid,
  content text,
  similarity float,
  metadata jsonb
)
language plpgsql
as $$
begin
  return query
  select
    g.id,
    e.chunk_content as content,
    1 - (e.embedding <=> query_embedding) as similarity,
    g.properties as metadata
  from golden_set_embeddings e
  join golden_set_library g on e.source_id = g.id
  where 1 - (e.embedding <=> query_embedding) > match_threshold
  and g.status = 'approved' -- 승인된 데이터만 검색!
  order by similarity desc
  limit match_count;
end;
$$;

-- 8. Row Level Security (RLS) - 보안 설정
alter table profiles enable row level security;
alter table projects enable row level security;
alter table golden_set_library enable row level security;

-- 내 프로필은 나만 본다
drop policy if exists "Users can view own profile" on profiles;
create policy "Users can view own profile" 
on profiles for select using (auth.uid() = id);

-- 내 프로젝트는 나만 본다
drop policy if exists "Users can view own projects" on projects;
create policy "Users can view own projects" 
on projects for select using (auth.uid() = user_id);

-- 승인된 골든셋은 누구나 본다 (AI Agent용)
drop policy if exists "Public read approved golden set" on golden_set_library;
create policy "Public read approved golden set"
on golden_set_library for select
using (status = 'approved');

-- ==========================================
-- 2026-01-19 Admin Page Overhaul Updates
-- ==========================================

-- 9. Knowledge Base (비정형 데이터)
create table if not exists public.knowledge_base (
  id bigint primary key generated always as identity,
  source_type text,
  title text,
  summary text,
  content text,
  relevance_score float,
  source_tier int,
  ai_reasoning text,
  rag_status text default 'pending',
  created_at timestamptz default now()
);

-- 10. Knowledge Base Embeddings
create table if not exists public.kb_embeddings (
  id bigint primary key generated always as identity,
  kb_id bigint references public.knowledge_base(id) on delete cascade,
  embedding vector(1536)
);

-- 11. Agent Prompts (AI Tuning)
create table if not exists public.agent_prompts (
  id uuid default gen_random_uuid() primary key,
  agent_id text not null,
  version text not null,
  system_prompt text not null,
  is_active boolean default false,
  is_live boolean default false,
  created_at timestamptz default now(),
  updated_at timestamptz default now()
);

-- 12. Schema Updates for Existing Tables
-- (Note: In a fresh init, these could be part of the CREATE TABLE, but keeping as ALTER for compatibility)
do $$
begin
  -- Golden Set Library Columns
  if not exists (select 1 from information_schema.columns where table_name = 'golden_set_library' and column_name = 'outcome_type') then
    alter table public.golden_set_library add column outcome_type text check (outcome_type in ('Success', 'Failure', 'Terminated'));
  end if;
  if not exists (select 1 from information_schema.columns where table_name = 'golden_set_library' and column_name = 'failure_reason') then
    alter table public.golden_set_library add column failure_reason text;
  end if;
  if not exists (select 1 from information_schema.columns where table_name = 'golden_set_library' and column_name = 'ip_status') then
    alter table public.golden_set_library add column ip_status text default 'Unknown';
  end if;
  if not exists (select 1 from information_schema.columns where table_name = 'golden_set_library' and column_name = 'smiles_code') then
    alter table public.golden_set_library add column smiles_code text;
  end if;
  if not exists (select 1 from information_schema.columns where table_name = 'golden_set_library' and column_name = 'patent_id') then
    alter table public.golden_set_library add column patent_id text;
  end if;
  if not exists (select 1 from information_schema.columns where table_name = 'golden_set_library' and column_name = 'patent_expiry') then
    alter table public.golden_set_library add column patent_expiry date;
  end if;

  -- Projects Columns
  if not exists (select 1 from information_schema.columns where table_name = 'projects' and column_name = 'error_detail') then
    alter table public.projects add column error_detail text;
  end if;
end $$;

-- 13. RLS for New Tables
alter table knowledge_base enable row level security;
alter table agent_prompts enable row level security;

create policy "Public read indexed knowledge base" on knowledge_base for select using (rag_status = 'indexed');
create policy "Public read active prompts" on agent_prompts for select using (is_active = true or is_live = true);

