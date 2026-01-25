-- Create commercial_reagents table if not exists
create table if not exists public.commercial_reagents (
  id uuid not null default gen_random_uuid (),
  ambeed_cat_no text not null,
  product_name text null,
  category text null,
  cas_number text null,
  smiles_code text null,
  molecular_weight text null,
  formula text null,
  price_data jsonb null,
  stock_status text null,
  product_url text null,
  crawled_at timestamp with time zone null default now(),
  target text null,
  properties jsonb null,
  source_name text null,
  summary text null,
  embedding public.vector null,
  ai_refined boolean null default false,
  is_manual_override boolean null default false,
  binding_affinity text null,
  isotype text null,
  host_species text null,
  orr_pct text null,
  os_months text null,
  pfs_months text null,
  constraint commercial_reagents_pkey primary key (id),
  constraint commercial_reagents_ambeed_cat_no_key unique (ambeed_cat_no)
) TABLESPACE pg_default;

-- Create indexes
create index IF not exists idx_commercial_reagents_category on public.commercial_reagents using btree (category) TABLESPACE pg_default;
create index IF not exists idx_commercial_reagents_cas on public.commercial_reagents using btree (cas_number) TABLESPACE pg_default;
create index IF not exists idx_commercial_reagents_source_name on public.commercial_reagents using btree (source_name) TABLESPACE pg_default;
create index IF not exists idx_commercial_reagents_embedding on public.commercial_reagents using hnsw (embedding vector_cosine_ops) TABLESPACE pg_default;
create index IF not exists idx_commercial_reagents_ai_refined on public.commercial_reagents using btree (ai_refined) TABLESPACE pg_default;
create index IF not exists idx_commercial_reagents_target on public.commercial_reagents using btree (target) TABLESPACE pg_default;
