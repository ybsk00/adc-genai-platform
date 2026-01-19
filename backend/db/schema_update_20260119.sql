-- 2026-01-19 Admin Page Overhaul Schema Update (Fixed)

-- 1. Knowledge Base Table (비정형 데이터)
create table if not exists public.knowledge_base (
  id bigint primary key generated always as identity,
  source_type text, -- 'News', 'PubMed', 'Report'
  title text,
  summary text,
  content text, -- 원문
  relevance_score float, -- AI 중요도 점수
  source_tier int, -- 1(High) ~ 3(Low)
  ai_reasoning text, -- 추천 사유
  rag_status text default 'pending', -- pending, indexed, excluded
  created_at timestamptz default now()
);

-- 2. Knowledge Base Embeddings (Vector DB)
create table if not exists public.kb_embeddings (
  id bigint primary key generated always as identity,
  kb_id bigint references public.knowledge_base(id) on delete cascade,
  embedding vector(1536) -- OpenAI text-embedding-3-small
);

-- 3. Golden Set Library Updates
-- 기존 테이블에 컬럼 추가 (DO 블록으로 안전하게 처리)
do $$
begin
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
end $$;

-- 4. Projects Table Update
do $$
begin
  if not exists (select 1 from information_schema.columns where table_name = 'projects' and column_name = 'error_detail') then
    alter table public.projects add column error_detail text;
  end if;
end $$;

-- 5. Agent Prompts Table (AI Tuning)
create table if not exists public.agent_prompts (
  id uuid default gen_random_uuid() primary key,
  agent_id text not null, -- 'structure', 'toxicology', etc.
  version text not null, -- 'v1.0', 'v1.1'
  system_prompt text not null,
  created_at timestamptz default now(),
  updated_at timestamptz default now()
);

-- Safely add columns to agent_prompts if they don't exist
do $$
begin
  if not exists (select 1 from information_schema.columns where table_name = 'agent_prompts' and column_name = 'is_active') then
    alter table public.agent_prompts add column is_active boolean default false;
  end if;
  if not exists (select 1 from information_schema.columns where table_name = 'agent_prompts' and column_name = 'is_live') then
    alter table public.agent_prompts add column is_live boolean default false;
  end if;
end $$;

-- 6. RLS Policies for new tables
alter table knowledge_base enable row level security;
alter table agent_prompts enable row level security;

-- 정책이 이미 존재할 수 있으므로 삭제 후 재생성 (또는 IF NOT EXISTS 체크)
drop policy if exists "Public read indexed knowledge base" on knowledge_base;
create policy "Public read indexed knowledge base"
on knowledge_base for select
using (rag_status = 'indexed');

drop policy if exists "Public read active prompts" on agent_prompts;
create policy "Public read active prompts"
on agent_prompts for select
using (is_active = true or is_live = true);
