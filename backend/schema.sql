-- Create golden_set_library table
create table if not exists public.golden_set_library (
  id uuid default gen_random_uuid() primary key,
  name text not null,
  category text,
  description text,
  properties jsonb default '{}'::jsonb,
  status text default 'draft', -- 'draft', 'approved', 'rejected'
  enrichment_source text,
  raw_data jsonb default '{}'::jsonb,
  outcome_type text, -- 'Success', 'Failure', 'Ongoing', 'Unknown'
  failure_reason text,
  created_at timestamp with time zone default timezone('utc'::text, now()) not null,
  updated_at timestamp with time zone default timezone('utc'::text, now()) not null
);

-- Enable Row Level Security (RLS)
alter table public.golden_set_library enable row level security;

-- Create policies
-- Allow read access to authenticated users
create policy "Allow read access for authenticated users"
on public.golden_set_library
for select
to authenticated
using (true);

-- Allow insert/update/delete access to service_role (backend) and admins
-- Note: service_role bypasses RLS by default, but explicit policies help if using authenticated client
create policy "Allow all access for service_role"
on public.golden_set_library
for all
to service_role
using (true)
with check (true);

-- Create index for faster lookups
create index if not exists idx_golden_set_library_status on public.golden_set_library(status);
create index if not exists idx_golden_set_library_name on public.golden_set_library(name);
