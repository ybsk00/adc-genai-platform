-- Fix Schema Relationships for Admin Dashboard
-- Run this in Supabase SQL Editor

-- 1. Ensure profiles table exists (Mirror of auth.users)
create table if not exists public.profiles (
  id uuid references auth.users on delete cascade primary key,
  email text,
  full_name text,
  credits int default 0,
  created_at timestamptz default now(),
  updated_at timestamptz default now()
);

-- 2. Trigger to automatically create profile on signup
create or replace function public.handle_new_user()
returns trigger as $$
begin
  insert into public.profiles (id, email, full_name)
  values (new.id, new.email, new.raw_user_meta_data->>'full_name');
  return new;
end;
$$ language plpgsql security definer;

-- Drop trigger if exists to avoid duplication error on create
drop trigger if exists on_auth_user_created on auth.users;
create trigger on_auth_user_created
  after insert on auth.users
  for each row execute procedure public.handle_new_user();

-- 3. Ensure projects table exists and has FK to profiles
create table if not exists public.projects (
  id uuid default gen_random_uuid() primary key,
  user_id uuid references public.profiles(id) on delete cascade not null,
  title text,
  status text default 'pending', -- pending, processing, completed, failed
  input_params jsonb,
  result_summary jsonb,
  report_url text,
  error_detail text,
  created_at timestamptz default now(),
  updated_at timestamptz default now()
);

-- 4. Add Foreign Key if it doesn't exist (Critical for join queries)
do $$
begin
  if not exists (
    select 1 from information_schema.table_constraints 
    where constraint_name = 'projects_user_id_fkey'
  ) then
    alter table public.projects 
    add constraint projects_user_id_fkey 
    foreign key (user_id) 
    references public.profiles(id) 
    on delete cascade;
  end if;
end $$;

-- 5. RLS Policies (Allow Admin to view all)
alter table public.profiles enable row level security;
alter table public.projects enable row level security;

-- Profiles: Users can view own, Admin can view all
drop policy if exists "Users can view own profile" on public.profiles;
create policy "Users can view own profile" on public.profiles
  for select using (auth.uid() = id);

-- Projects: Users can view own, Admin can view all
drop policy if exists "Users can view own projects" on public.projects;
create policy "Users can view own projects" on public.projects
  for select using (auth.uid() = user_id);

-- Note: Service Role (Backend) bypasses RLS, so these are for Frontend direct access if needed.
