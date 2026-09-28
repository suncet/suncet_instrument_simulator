-- Run this in Supabase's SQL Editor, then run catalog.sql. Both are repeatable.
begin;

create table if not exists public.colortable_studies (
    id text primary key,
    title text not null,
    is_open boolean not null default true
);
create table if not exists public.colortable_palettes (
    study_id text not null references public.colortable_studies(id),
    slug text not null,
    number integer not null,
    title text not null,
    primary key (study_id, slug),
    unique (study_id, number)
);
create table if not exists public.colortable_favorites (
    study_id text not null,
    palette_slug text not null,
    user_id uuid not null references auth.users(id) on delete cascade,
    created_at timestamptz not null default now(),
    primary key (study_id, user_id, palette_slug),
    foreign key (study_id, palette_slug) references public.colortable_palettes(study_id, slug)
);
create index if not exists colortable_favorites_counts
    on public.colortable_favorites(study_id, palette_slug);
create table if not exists public.colortable_vote_limits (
    user_id uuid primary key references auth.users(id) on delete cascade,
    window_start timestamptz not null,
    requests integer not null
);

alter table public.colortable_studies enable row level security;
alter table public.colortable_palettes enable row level security;
alter table public.colortable_favorites enable row level security;
alter table public.colortable_vote_limits enable row level security;

revoke all on public.colortable_studies, public.colortable_palettes,
    public.colortable_favorites, public.colortable_vote_limits from public, anon, authenticated;
grant select on public.colortable_studies, public.colortable_palettes to anon, authenticated;
grant select on public.colortable_favorites to authenticated;

drop policy if exists colortable_read_studies on public.colortable_studies;
create policy colortable_read_studies on public.colortable_studies for select to anon, authenticated using (true);
drop policy if exists colortable_read_palettes on public.colortable_palettes;
create policy colortable_read_palettes on public.colortable_palettes for select to anon, authenticated using (true);
drop policy if exists colortable_read_own on public.colortable_favorites;
create policy colortable_read_own on public.colortable_favorites for select to authenticated
    using (user_id = (select auth.uid()));

-- Only this function can mutate votes. Caller identity always comes from Auth.
create or replace function public.set_colortable_favorite(p_study text, p_palette text, p_favorite boolean)
returns boolean
language plpgsql
security definer
set search_path = ''
as $$
declare
    voter uuid := auth.uid();
    request_count integer;
begin
    if voter is null then
        raise exception 'Sign in anonymously before voting.' using errcode = '42501';
    end if;
    if p_favorite is null then
        raise exception 'A favorite state is required.' using errcode = '22023';
    end if;
    -- The shared lock also serializes a concurrent owner request to close voting.
    perform 1 from public.colortable_studies where id = p_study and is_open for share;
    if not found then
        raise exception 'Voting is closed or the study does not exist.' using errcode = '22023';
    end if;
    perform 1 from public.colortable_palettes where study_id = p_study and slug = p_palette;
    if not found then
        raise exception 'Unknown palette.' using errcode = '22023';
    end if;
    insert into public.colortable_vote_limits as limits (user_id, window_start, requests)
        values (voter, now(), 1)
        on conflict (user_id) do update set
            window_start = case when limits.window_start < now() - interval '1 minute' then now() else limits.window_start end,
            requests = case when limits.window_start < now() - interval '1 minute' then 1 else limits.requests + 1 end
        returning requests into request_count;
    if request_count > 90 then
        raise exception 'Too many vote changes. Please wait a minute.' using errcode = 'P0001';
    end if;
    if p_favorite then
        insert into public.colortable_favorites(study_id, palette_slug, user_id)
            values (p_study, p_palette, voter) on conflict do nothing;
    else
        delete from public.colortable_favorites
            where study_id = p_study and palette_slug = p_palette and user_id = voter;
    end if;
    return p_favorite;
end;
$$;

-- Aggregate counts are public; individual voter IDs and selections are not.
create or replace function public.colortable_counts(p_study text)
returns table(palette_slug text, favorites bigint)
language sql
stable
security definer
set search_path = ''
as $$
    select p.slug, count(f.user_id)
    from public.colortable_palettes p
    left join public.colortable_favorites f on f.study_id = p.study_id and f.palette_slug = p.slug
    where p.study_id = p_study
    group by p.slug;
$$;

revoke all on function public.set_colortable_favorite(text, text, boolean) from public, anon, authenticated;
grant execute on function public.set_colortable_favorite(text, text, boolean) to authenticated;
revoke all on function public.colortable_counts(text) from public, anon, authenticated;
grant execute on function public.colortable_counts(text) to anon, authenticated;

commit;
