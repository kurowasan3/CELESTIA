-- CELESTIA：学習記録の同期用（Supabase の SQL Editor で実行する。何度実行してもよい）
-- Supabase プロジェクトは Lemma と共通（アカウントも共通）。記録は Lemma（lemma_state）とは別の表に1人1行。
-- 各ユーザーは自分の行だけを読み書きできる（行単位のアクセス制限 RLS）。
-- ここで作る・変える名前（表・制約・ポリシー）はすべて CELESTIA のもので、Lemma の表・関数・ポリシーには触れない。
-- このあとに supabase/guard.sql（サーバー側の書き込み制限）を実行する（docs/data-integrity.md §6.3）。
create table if not exists public.celestia_state (
  user_id    uuid primary key references auth.users(id) on delete cascade,
  state      jsonb not null,
  rev        integer not null default 1,
  updated_at timestamptz not null default now(),
  constraint celestia_state_size check (pg_column_size(state) < 5000000)
);

-- v0.5 までの schema.sql で作った表は、大きさの制約の名前が state_size（Lemma と同じ名前）だったので、CELESTIA の名前に変える（中身は同じ）
do $$ begin
  if exists (select 1 from pg_constraint where conname = 'state_size' and conrelid = 'public.celestia_state'::regclass) then
    alter table public.celestia_state rename constraint state_size to celestia_state_size;
  end if;
end $$;

alter table public.celestia_state enable row level security;

-- ポリシーも CELESTIA の名前にする（v0.5 までの名前は Lemma と同じだったので、あれば消してから作り直す。celestia_state のポリシーだけ）
drop policy if exists "自分の記録を読む" on public.celestia_state;
drop policy if exists "自分の記録を作る" on public.celestia_state;
drop policy if exists "自分の記録を更新" on public.celestia_state;
drop policy if exists "自分の記録を消す" on public.celestia_state;
drop policy if exists "celestia_state: 自分の記録を読む" on public.celestia_state;
drop policy if exists "celestia_state: 自分の記録を作る" on public.celestia_state;
drop policy if exists "celestia_state: 自分の記録を更新" on public.celestia_state;
drop policy if exists "celestia_state: 自分の記録を消す" on public.celestia_state;
create policy "celestia_state: 自分の記録を読む" on public.celestia_state
  for select using (auth.uid() = user_id);
create policy "celestia_state: 自分の記録を作る" on public.celestia_state
  for insert with check (auth.uid() = user_id);
create policy "celestia_state: 自分の記録を更新" on public.celestia_state
  for update using (auth.uid() = user_id) with check (auth.uid() = user_id);
create policy "celestia_state: 自分の記録を消す" on public.celestia_state
  for delete using (auth.uid() = user_id);

-- 退会用の関数 public.delete_my_account() は Lemma の schema.sql ですでに作成済み（共通で使う。ここでは作らない・変えない）。
-- アカウントごと消えるので、Lemma の記録（lemma_state）も一緒に消える（celestia_state は on delete cascade で消える）。
