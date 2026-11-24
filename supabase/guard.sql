-- CELESTIA v0.6：学習記録をサーバー側でも守る（schema.sql のあとに、Supabase の SQL Editor で実行する）
--   Lemma v0.26 の supabase/guard.sql を、CELESTIA の表 celestia_state 用にしたもの。
--   ・Supabase プロジェクトは Lemma と共通。ここで作る schema・関数・トリガーの名前はすべて celestia_ で始まり、
--     Lemma の表 lemma_state・schema lemma_private・トリガー lemma_state_guard には一切触れない（tests/server_guard.test.mjs で確認）。
--   ・何度実行してもよい（関数は create or replace、トリガーは作り直すだけ）。表と既存の記録は変えない。
--   ・適用前後の確認・戻し方は docs/data-integrity.md §6.3。
--
-- アプリ（クライアント）の合わせ方（js/store.js mergeState）は「和集合・新しい方・多い方」なので、通常の同期で
-- 学習記録が減ったり古くなったりすることは無い。減ってよいのは「消した印」resetAt を進めた明示的な消去・置き換えだけ。
-- そこで、表 celestia_state への書き込みごとに、次の規則を確かめる（破れば SQLSTATE 'CE0xx' のエラーで断り、表は変えない）。
--
--   CE001 形：state はオブジェクト、log は配列で各要素は {p: 文字列, t: 数値}、pstat などはオブジェクト、predictions は配列で各要素は {p, t}
--   CE002 版番号：作成は rev = 1、更新は rev = 前の rev + 1
--   CE003 記録の版：state.schema は 1〜100 の整数で、更新で下げない（v0.5 までのアプリが schema 1 で書き戻すのも断る）
--   CE004 解答の記録：前の記録（問題ID p と時刻 t の組）がすべて残っている（resetAt を進めた消去なら、その時刻より前は消えてよい）
--   CE005 解答の記録の改変：残っている記録から項目を消していない（値の更新＝再採点・メモ・異議は許す）
--   CE006 消した印 resetAt：数値で、戻さない。進めるときはサーバーの時計で1日より先の未来にしない
--   CE007 学習状態：pstat・problems・bookmarks・explains・misc の項目が消えていない・時刻（と problems・misc の回数）が戻っていない
--         （resetAt を進めた消去なら、その時刻より前のものは消えてよい）
--   CE008 重複：同じ (p, t) の解答の記録を増やしていない
--   CE009 同じ時刻のままの書き換え：
--         pstat：値が変わるなら、新しい値は正しい範囲（S>0、1≦D≦10、reps・lapses は 0 以上の整数で lapses≦reps）。
--                last が同じなら reps を減らさない。reps も同じなら lapses の差が1まで（再採点）
--         bookmarks：t が同じなら値を変えない（付け外しは新しい t になる）
--         explains：前にあったメモ（t と本文）は、新しい5件から押し出された場合を除いて残る（同じ t のまま本文を変えない）
--   CE010 予測の記録：前の予測（p と t の組）がすべて残っている（resetAt を進めた消去なら、その時刻より前は消えてよい）
--   localResetAt（端末だけの印）は、送られてきても保存しない（取り除く）。
--
-- 守らないもの（通常の操作で正しく変わるため）：settings・practice・simPrefs（最後に変えた方）、streak（日をあけると戻る）、
--   explains の古いメモ（問題ごとに新しい5件だけ残す）、misc の状態（cleared など。新しい方に入れ替わる）。

-- 関数は API に出ない schema（celestia_private）に置く（public の関数は Supabase の API から呼べてしまうため）。
-- トリガーは書き込んだ利用者の権限で動くので、authenticated にこの schema の利用を許す。
create schema if not exists celestia_private;
revoke all on schema celestia_private from public, anon;
grant usage on schema celestia_private to authenticated, service_role;

-- 入っている guard.sql の版（適用済みかの確認用：select celestia_private.guard_version();）
create or replace function celestia_private.guard_version() returns text
language sql immutable set search_path = '' as $$ select '0.6.0'::text $$;

create or replace function celestia_private.celestia_num(v jsonb) returns numeric
language sql immutable set search_path = '' as $$
  select case when jsonb_typeof(v) = 'number' then (v #>> '{}')::numeric end
$$;

-- 対応表（問題ID → {時刻 field, …}）で、守る項目（時刻 >= cutoff）が消えた・時刻が戻ったものの数と例
create or replace function celestia_private.celestia_map_lost(o jsonb, n jsonb, field text, cutoff numeric, out lost bigint, out example text)
language sql immutable set search_path = '' as $$
  select count(*), min(a.key)
  from jsonb_each(case when jsonb_typeof(o) = 'object' then o else '{}'::jsonb end) a
  left join jsonb_each(case when jsonb_typeof(n) = 'object' then n else '{}'::jsonb end) b on b.key = a.key
  where coalesce(celestia_private.celestia_num(a.value -> field), 0) >= cutoff
    and (b.key is null
         or (celestia_private.celestia_num(a.value -> field) is not null
             and coalesce(celestia_private.celestia_num(b.value -> field), -1) < celestia_private.celestia_num(a.value -> field)))
$$;

-- 復習状態の値（{S, D, last, reps, lapses}）が正しい範囲か
create or replace function celestia_private.celestia_review_ok(v jsonb) returns boolean
language sql immutable set search_path = '' as $$
  select jsonb_typeof(v) = 'object'
    and coalesce(celestia_private.celestia_num(v -> 'S') > 0, false)
    and coalesce(celestia_private.celestia_num(v -> 'D') between 1 and 10, false)
    and celestia_private.celestia_num(v -> 'last') is not null
    and coalesce(celestia_private.celestia_num(v -> 'reps') >= 0 and celestia_private.celestia_num(v -> 'reps') = floor(celestia_private.celestia_num(v -> 'reps')), false)
    and coalesce(celestia_private.celestia_num(v -> 'lapses') >= 0 and celestia_private.celestia_num(v -> 'lapses') = floor(celestia_private.celestia_num(v -> 'lapses'))
                 and celestia_private.celestia_num(v -> 'lapses') <= celestia_private.celestia_num(v -> 'reps'), false)
$$;

-- 復習状態の対応表で、前からある項目（last >= cutoff）の値が変わったもののうち、不正なものの数と例。
--   regrade = true（pstat）：last と reps が同じなら lapses の差1までの変化を許す（再採点は前の状態から計算し直すため S・D も変わる）
--   regrade = false：last と reps が同じなら値を変えない（Lemma の解法カード用。CELESTIA では今は使わない）
--   last が進んだ変化は、値が正しい範囲なら許す（復習・2台の復習の計算し直し。古い版のアプリは reps が減る合わせ方をすることがある）
create or replace function celestia_private.celestia_review_bad(o jsonb, n jsonb, cutoff numeric, regrade boolean, out bad bigint, out example text)
language sql immutable set search_path = '' as $$
  select count(*), min(a.key)
  from jsonb_each(case when jsonb_typeof(o) = 'object' then o else '{}'::jsonb end) a
  join jsonb_each(case when jsonb_typeof(n) = 'object' then n else '{}'::jsonb end) b on b.key = a.key
  where coalesce(celestia_private.celestia_num(a.value -> 'last'), 0) >= cutoff
    and b.value is distinct from a.value
    and (not celestia_private.celestia_review_ok(b.value)
         or (celestia_private.celestia_num(a.value -> 'last') = celestia_private.celestia_num(b.value -> 'last')
             and celestia_private.celestia_num(a.value -> 'reps') is not null
             and (celestia_private.celestia_num(b.value -> 'reps') < celestia_private.celestia_num(a.value -> 'reps')
                  or (celestia_private.celestia_num(b.value -> 'reps') = celestia_private.celestia_num(a.value -> 'reps')
                      and (not regrade
                           or abs(celestia_private.celestia_num(b.value -> 'lapses') - coalesce(celestia_private.celestia_num(a.value -> 'lapses'), 0)) > 1)))))
$$;

create or replace function celestia_private.celestia_state_guard() returns trigger
language plpgsql set search_path = '' as $$
declare
  o jsonb; n jsonb;
  old_reset numeric; new_reset numeric; cutoff numeric;
  old_schema numeric; new_schema numeric;
  cnt bigint; ex text; m text; old_n bigint;
begin
  -- 端末だけの印は保存しない
  if jsonb_typeof(new.state) = 'object' and new.state ? 'localResetAt' then
    new.state := new.state - 'localResetAt';
  end if;
  n := new.state;
  ---------------------------------------------------------------- 形（作成・更新とも）
  if jsonb_typeof(n) is distinct from 'object' or jsonb_typeof(n -> 'log') is distinct from 'array' then
    raise exception 'celestia: state はオブジェクト、state.log は配列でなければなりません' using errcode = 'CE001';
  end if;
  select count(*) into cnt from jsonb_array_elements(n -> 'log') e
    where jsonb_typeof(e) <> 'object' or jsonb_typeof(e -> 'p') is distinct from 'string' or jsonb_typeof(e -> 't') is distinct from 'number';
  if cnt > 0 then
    raise exception 'celestia: 問題ID p（文字列）と時刻 t（数値）の無い解答の記録が % 件あります', cnt using errcode = 'CE001';
  end if;
  foreach m in array array['pstat', 'problems', 'bookmarks', 'explains', 'misc'] loop
    if n ? m and jsonb_typeof(n -> m) is distinct from 'object' then
      raise exception 'celestia: state.% はオブジェクトでなければなりません', m using errcode = 'CE001';
    end if;
  end loop;
  if n ? 'predictions' then
    if jsonb_typeof(n -> 'predictions') is distinct from 'array' then
      raise exception 'celestia: state.predictions は配列でなければなりません' using errcode = 'CE001';
    end if;
    select count(*) into cnt from jsonb_array_elements(n -> 'predictions') e
      where jsonb_typeof(e) <> 'object' or jsonb_typeof(e -> 'p') is distinct from 'string' or jsonb_typeof(e -> 't') is distinct from 'number';
    if cnt > 0 then
      raise exception 'celestia: 問題ID p（文字列）と時刻 t（数値）の無い予測の記録が % 件あります', cnt using errcode = 'CE001';
    end if;
  end if;
  new_schema := celestia_private.celestia_num(n -> 'schema');
  if new_schema is null or new_schema <> floor(new_schema) or new_schema < 1 or new_schema > 100 then
    raise exception 'celestia: state.schema は 1〜100 の整数でなければなりません（%）', coalesce(n ->> 'schema', '無し') using errcode = 'CE003';
  end if;
  if n ? 'resetAt' and celestia_private.celestia_num(n -> 'resetAt') is null then
    raise exception 'celestia: resetAt は数値でなければなりません' using errcode = 'CE006';
  end if;
  new_reset := coalesce(celestia_private.celestia_num(n -> 'resetAt'), 0);

  if tg_op = 'INSERT' then
    if new.rev is distinct from 1 then
      raise exception 'celestia: 作成時の rev は 1 でなければなりません（%）', new.rev using errcode = 'CE002';
    end if;
    if new_reset > extract(epoch from now()) * 1000 + 86400000 then
      raise exception 'celestia: resetAt が未来です' using errcode = 'CE006';
    end if;
    select count(*) - count(distinct (e ->> 'p', celestia_private.celestia_num(e -> 't'))) into cnt from jsonb_array_elements(n -> 'log') e;
    if cnt > 0 then
      raise exception 'celestia: 同じ問題・同じ時刻の解答の記録が % 件重複しています', cnt using errcode = 'CE008';
    end if;
    return new;
  end if;

  ---------------------------------------------------------------- 更新
  o := old.state;
  if new.rev is distinct from old.rev + 1 then
    raise exception 'celestia: rev は % のはずです（送られたのは %）', old.rev + 1, new.rev using errcode = 'CE002';
  end if;
  old_schema := coalesce(celestia_private.celestia_num(o -> 'schema'), 1);
  if new_schema < old_schema then
    raise exception 'celestia: 記録の版を % から % に下げる書き込みです', old_schema, new_schema using errcode = 'CE003';
  end if;
  old_reset := coalesce(celestia_private.celestia_num(o -> 'resetAt'), 0);
  if new_reset < old_reset then
    raise exception 'celestia: 消した印 resetAt を戻す書き込みです' using errcode = 'CE006';
  end if;
  if new_reset > old_reset then
    if new_reset > extract(epoch from now()) * 1000 + 86400000 then
      raise exception 'celestia: resetAt が未来です（端末の時計を確かめてください）' using errcode = 'CE006';
    end if;
    cutoff := new_reset;        -- 明示的な消去・置き換え：この時刻より前の記録は消えてよい
  else
    cutoff := '-infinity';      -- 通常の同期：何も消えてはいけない
  end if;

  -- 解答の記録：前の記録（形の正しいもの）が残っているか・項目が消えていないか・重複を増やしていないか
  with a as (
    select e ->> 'p' as p, celestia_private.celestia_num(e -> 't') as t, e from jsonb_array_elements(
      case when jsonb_typeof(o -> 'log') = 'array' then o -> 'log' else '[]'::jsonb end) e
    where jsonb_typeof(e) = 'object' and jsonb_typeof(e -> 'p') = 'string' and jsonb_typeof(e -> 't') = 'number'
  ), b as (
    select distinct on (e ->> 'p', celestia_private.celestia_num(e -> 't')) e ->> 'p' as p, celestia_private.celestia_num(e -> 't') as t, e
    from jsonb_array_elements(n -> 'log') e
  )
  select count(*), min(a.p || '@' || a.t) into cnt, ex from a left join b on b.p = a.p and b.t = a.t
    where a.t >= cutoff and b.p is null;
  if cnt > 0 then
    raise exception 'celestia: 前にあった解答の記録 % 件が無くなる書き込みです（例 %）', cnt, ex using errcode = 'CE004';
  end if;
  with a as (
    select e ->> 'p' as p, celestia_private.celestia_num(e -> 't') as t, e from jsonb_array_elements(
      case when jsonb_typeof(o -> 'log') = 'array' then o -> 'log' else '[]'::jsonb end) e
    where jsonb_typeof(e) = 'object' and jsonb_typeof(e -> 'p') = 'string' and jsonb_typeof(e -> 't') = 'number'
  ), b as (
    select e ->> 'p' as p, celestia_private.celestia_num(e -> 't') as t, e from jsonb_array_elements(n -> 'log') e
  )
  select count(*), min(a.p || '@' || a.t) into cnt, ex from a join b on b.p = a.p and b.t = a.t
    where a.t >= cutoff and exists (select 1 from jsonb_object_keys(a.e) k where not (b.e ? k));   -- 消した印（resetAt）より前の記録は CE004 と同じく対象外
  if cnt > 0 then
    raise exception 'celestia: 解答の記録から項目を消す書き込みです（例 %）', ex using errcode = 'CE005';
  end if;
  select (select count(*) - count(distinct (e ->> 'p', celestia_private.celestia_num(e -> 't'))) from jsonb_array_elements(n -> 'log') e)
       - (select count(*) - count(distinct (e ->> 'p', celestia_private.celestia_num(e -> 't'))) from jsonb_array_elements(
            case when jsonb_typeof(o -> 'log') = 'array' then o -> 'log' else '[]'::jsonb end) e
          where jsonb_typeof(e) = 'object' and jsonb_typeof(e -> 'p') = 'string' and jsonb_typeof(e -> 't') = 'number')
    into cnt;
  if cnt > 0 then
    raise exception 'celestia: 同じ問題・同じ時刻の解答の記録を % 件増やす書き込みです', cnt using errcode = 'CE008';
  end if;

  -- 学習状態：消えていない・時刻が戻っていない（problems は回数も）
  foreach m in array array['pstat', 'problems', 'misc'] loop
    select l.lost, l.example into cnt, ex from celestia_private.celestia_map_lost(o -> m, n -> m, 'last', cutoff) l;
    if cnt > 0 then
      raise exception 'celestia: % の % 件が消える・古くなる書き込みです（例 %）', m, cnt, ex using errcode = 'CE007';
    end if;
  end loop;
  select l.lost, l.example into cnt, ex from celestia_private.celestia_map_lost(o -> 'bookmarks', n -> 'bookmarks', 't', cutoff) l;
  if cnt > 0 then
    raise exception 'celestia: bookmarks の % 件が消える・古くなる書き込みです（例 %。外すときは on:false と新しい時刻で残します）', cnt, ex using errcode = 'CE007';
  end if;
  select count(*), min(a.key) into cnt, ex
    from jsonb_each(case when jsonb_typeof(o -> 'problems') = 'object' then o -> 'problems' else '{}'::jsonb end) a
    join jsonb_each(n -> 'problems') b on b.key = a.key
    where coalesce(celestia_private.celestia_num(a.value -> 'last'), 0) >= cutoff
      and celestia_private.celestia_num(a.value -> 'n') is not null
      and coalesce(celestia_private.celestia_num(b.value -> 'n'), -1) < celestia_private.celestia_num(a.value -> 'n');
  if cnt > 0 then
    raise exception 'celestia: problems の解いた回数が減る書き込みです（例 %）', ex using errcode = 'CE007';
  end if;
  select count(*), min(a.key) into cnt, ex
    from jsonb_each(case when jsonb_typeof(o -> 'misc') = 'object' then o -> 'misc' else '{}'::jsonb end) a
    join jsonb_each(case when jsonb_typeof(n -> 'misc') = 'object' then n -> 'misc' else '{}'::jsonb end) b on b.key = a.key
    where coalesce(celestia_private.celestia_num(a.value -> 'last'), 0) >= cutoff
      and celestia_private.celestia_num(a.value -> 'n') is not null
      and coalesce(celestia_private.celestia_num(b.value -> 'n'), -1) < celestia_private.celestia_num(a.value -> 'n');
  if cnt > 0 then
    raise exception 'celestia: misc（誤概念）の検出回数が減る書き込みです（例 %）', ex using errcode = 'CE007';
  end if;
  -- 予測の記録：前の予測（p と t の組）が残っているか
  with a as (
    select e ->> 'p' as p, celestia_private.celestia_num(e -> 't') as t from jsonb_array_elements(
      case when jsonb_typeof(o -> 'predictions') = 'array' then o -> 'predictions' else '[]'::jsonb end) e
    where jsonb_typeof(e) = 'object' and jsonb_typeof(e -> 'p') = 'string' and jsonb_typeof(e -> 't') = 'number'
  ), b as (
    select distinct e ->> 'p' as p, celestia_private.celestia_num(e -> 't') as t from jsonb_array_elements(
      case when jsonb_typeof(n -> 'predictions') = 'array' then n -> 'predictions' else '[]'::jsonb end) e
  )
  select count(*), min(a.p || '@' || a.t) into cnt, ex from a left join b on b.p = a.p and b.t = a.t
    where a.t >= cutoff and b.p is null;
  if cnt > 0 then
    raise exception 'celestia: 前にあった予測の記録 % 件が無くなる書き込みです（例 %）', cnt, ex using errcode = 'CE010';
  end if;
  -- 同じ時刻のままの書き換え（CE009）
  select l.bad, l.example into cnt, ex from celestia_private.celestia_review_bad(o -> 'pstat', n -> 'pstat', cutoff, true) l;
  if cnt > 0 then
    raise exception 'celestia: pstat の % 件が、正しくない値か、同じ時刻のまま書き換えられています（例 %）', cnt, ex using errcode = 'CE009';
  end if;
  select count(*), min(a.key) into cnt, ex
    from jsonb_each(case when jsonb_typeof(o -> 'bookmarks') = 'object' then o -> 'bookmarks' else '{}'::jsonb end) a
    join jsonb_each(n -> 'bookmarks') b on b.key = a.key
    where coalesce(celestia_private.celestia_num(a.value -> 't'), 0) >= cutoff
      and celestia_private.celestia_num(a.value -> 't') = celestia_private.celestia_num(b.value -> 't')
      and b.value is distinct from a.value;
  if cnt > 0 then
    raise exception 'celestia: bookmarks の % 件が、同じ時刻のまま書き換えられています（例 %。付け外しは新しい時刻になります）', cnt, ex using errcode = 'CE009';
  end if;
  -- メモ：前にあったメモ（t と本文）は残る。ただし新しい5件から押し出された（残ったメモがどれもそれより新しく、5件ある）なら外れてよい。
  --   同じ時刻のメモが別の本文で入っていれば、押し出しに見えても書き換え（CE009）
  with a as (
    select k.key, mm.value as memo, celestia_private.celestia_num(mm.value -> 't') as t
    from jsonb_each(case when jsonb_typeof(o -> 'explains') = 'object' then o -> 'explains' else '{}'::jsonb end) k,
         jsonb_array_elements(case when jsonb_typeof(k.value) = 'array' then k.value else '[]'::jsonb end) mm
    where celestia_private.celestia_num(mm.value -> 't') >= cutoff
  ), b as (
    select k.key, case when jsonb_typeof(k.value) = 'array' then k.value else '[]'::jsonb end as arr
    from jsonb_each(n -> 'explains') k
  ), gone as (
    select a.key, a.t, a.memo, b.arr from a left join b on b.key = a.key
    where b.arr is null
       or not exists (select 1 from jsonb_array_elements(b.arr) x where x -> 't' = a.memo -> 't' and x -> 'text' is not distinct from a.memo -> 'text')
  )
  select count(*) filter (where gone.arr is null
                             or jsonb_array_length(gone.arr) < 5
                             or exists (select 1 from jsonb_array_elements(gone.arr) x where coalesce(celestia_private.celestia_num(x -> 't'), 0) < gone.t)),
         min(gone.key || '@' || gone.t),
         count(*) filter (where gone.arr is not null and exists (select 1 from jsonb_array_elements(gone.arr) x where x -> 't' = gone.memo -> 't'))
    into cnt, ex, old_n
    from gone;
  if old_n > 0 then
    raise exception 'celestia: メモが同じ時刻のまま書き換えられています（例 %）', ex using errcode = 'CE009';
  end if;
  if cnt > 0 then
    raise exception 'celestia: explains（メモ）の % 件が消える書き込みです（例 %）', cnt, ex using errcode = 'CE007';
  end if;
  return new;
end $$;

drop trigger if exists celestia_state_guard on public.celestia_state;
create trigger celestia_state_guard before insert or update on public.celestia_state
  for each row execute function celestia_private.celestia_state_guard();

-- 関数は authenticated（アプリの利用者）と service_role（管理用の鍵）だけが使う。anon・PUBLIC からは使えない
revoke all on all functions in schema celestia_private from public, anon;
grant execute on all functions in schema celestia_private to authenticated, service_role;
