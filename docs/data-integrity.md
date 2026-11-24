# 学習記録の保全（v0.6）

学習記録を失わないための仕組み。Lemma v0.28（docs/data-integrity.md、v0.23〜v0.26）の仕組みを CELESTIA に引き継いだもの。
実装は `js/store.js`（保存・移行・点検・バックアップ・合わせ方）、`js/sync.js`（通信）、`js/app.js`（画面からの呼び出し）。
テストは `tests/store.test.mjs`・`tests/store_edge.test.mjs`・`tests/sync.test.mjs`・`tests/server_guard.test.mjs`・`tests/e2e.mjs`。
サーバー側の書き込み制限は `supabase/guard.sql`（Supabase で実行が必要。§6）。

Lemma との違い：解法カード（`cards`）・計算トレーニング（`drill`）・グラフの式（`graph*Exprs`）は無い。
代わりに CELESTIA の `predictions`（予測 → 実験の記録）・`misc`（誤概念の検出）・`simPrefs`（シミュレーションの条件）がある。
保存キー・控え・バックアップの形式名・サーバーの関数とエラーの code は、すべて CELESTIA だけのもの（同じドメインの Lemma・ぶつりば、同じ Supabase プロジェクトの Lemma と混ざらない）。

## 1. 保存しているもの

| 保存先（localStorage のキー） | 中身 | 同期 | バックアップ |
|---|---|---|---|
| `celestia.v1` | 学習記録の本体（下の表） | する | 含む |
| `celestia.snap.v1` | 自動の控え（最大2件。移行・読み込み・消去の前の `celestia.v1` の文字列） | しない | 含まない |
| `celestia.v1.broken.<時刻>` | 読めなかった記録の退避（元の文字列そのまま） | しない | 含まない（設定から個別に書き出せる） |
| `celestia.session.v1` | 中断した演習（14日で捨てる） | しない | 含まない |
| `celestia.auth` | ログイン情報（トークン） | — | 含まない（安全のため） |

IndexedDB・Cookie は使っていない。サーバー（Supabase。Lemma と共通のプロジェクト）は表 `celestia_state` に1人1行（`state` に `celestia.v1` と同じ形、`rev` に版番号）。

`celestia.v1` の中身（schema 5）：

| キー | 形 | 意味 |
|---|---|---|
| `log` | `[{ p, t, g, ok?, sg?, c?, h?, min?, in?, md?, note?, disp?, feel? }]` | 解答の記録。`p` 問題ID・`t` 時刻(ms)・`g` 評価1〜4・`md` 答えた形式（input・choice・self）。**`(p, t)` が1回の解答の識別子** |
| `pstat` | `{ 問題ID: { S, D, last, reps, lapses } }` | 問題ごとの復習状態（`scheduler.js` `review()`） |
| `problems` | `{ 問題ID: { last, n } }` | 最後に解いた時刻と回数 |
| `bookmarks` | `{ 問題ID: { on, t } }` | ブックマーク（外すと `on:false` を残す） |
| `explains` | `{ 問題ID: [{ t, text }] }` | 気づきのメモ（問題ごとに新しい5件） |
| `predictions` | `[{ p, t, predicted, answer?, conf? }]` | 予測 → 実験の記録（`(p, t)` が識別子） |
| `misc` | `{ 誤概念ID: { n, last, cleared? } }` | 誤概念ごとの検出回数と最後の時刻 |
| `streak` `csvExportedAt` | | 連続日数・CSV書き出し時刻 |
| `settings` `practice` `simPrefs` `settingsAt` | | 設定・演習の選択・シミュレーションの条件（最後に変えた方を採る） |
| `resetAt` | 時刻 | 「すべて消す」「置き換え」をした印（同期で古い記録を持ち込まない） |
| `localResetAt` | 時刻 | この端末で「ログアウトして消す」等をした印（同期しない。ほかのタブの古い記録を持ち込まない） |
| `schema` | | 記録の形の版 |

`predictions`・`misc` は v0.6 のアプリはまだ書き込まない（予測 → 実験の画面は未実装）。形・合わせ方・バックアップ・サーバーの保護は先に用意してある。

## 2. ローカル保存（`createStore().save()`）

Lemma と同じ。読み書きは「成功・無い・失敗」を区別し（`getItem` の例外を「記録が無い」とは扱わない）、

1. 保存を止めている理由（`blocked`）があれば書かない（`read` だけは読み直して解除を試す）。
2. 新しい schema（`newer`）・構造が正しくない（`invalid`）・JSON にできない（`serialize`）記録は書かない。
3. 保存されている記録がこのタブの最後に読み書きした文字列と違えば、壊れていれば退避、新しい schema なら止め、正しい記録なら `mergeTabs` で合わせる（同じなら書かない）。
4. `setItem`。容量不足なら古い控えから捨てて再試行、それでも入らなければ `quota`。どちらも前の値のまま。
5. 読み返して確かめる（別のタブが直後に書いたら合わせ直す。最大4回、収まらなければ `conflict`）。

失敗の理由は画面に1回知らせる（`app.js` `SAVE_MSG`）。保存している記録が300万文字を超えたら1回知らせる。
初めて解いたときに `navigator.storage.persist()` を頼む（Safari などが保存領域を勝手に消さないように。対応していない環境では何もしない）。

## 3. 読み込み（起動時の `load()`）

| 状態（status） | 処理 | 保存 |
|---|---|---|
| `new` 何も無い | 新しい記録 | する |
| `unavailable` 読めない（`getItem` の例外） | 新しい記録で画面を動かすが、記録が無いとは扱わない | 止める。読めたら保存されている記録と合わせてから保存 |
| `broken` JSON として読めない | 元の文字列を `celestia.v1.broken.<時刻>` に退避 | 退避できたらする |
| `salvaged` 一部の形が壊れている | 元を退避し、壊れた入れ物だけ初期値、ほかは使う | 退避できたらする |
| `migrated` 古い schema（v0.5 までの記録） | 移行（§4）。移行の前に控えを残す | する（次に保存するとき） |
| `migrate_failed` 移行の途中で例外 | 元を退避し、移行前の形のまま画面で使う | 止める（次の起動で再試行） |
| `newer` 新しい schema | 形を変えずに画面で使う | 止める |

いずれも起動時にバナーで知らせる（`migrated` は知らせない。移行の前の記録は設定の「自動の控え」から戻せる）。

## 4. schema と移行（v0.5 → v0.6）

- 今の版：`SCHEMA = 5`。読める最も古い版：`MIN_SCHEMA = 1`。schema の無い記録は 1 とみなす。
- **v0.5 までの記録は schema 1**（v0.1 から変わっていない）。ただし v0.5 の「読み込む」は、Lemma から写した名残で **schema 4** を書いていた。どちらも同じ形なので、1→2→3→4 は版を進めるだけ、4→5 が v0.6 への移行：
  - `practice.fields` の「分野横断」（`SYN-01`。v0.4 で分かれた）を `MSY-01`・`ESY-01` に読み替える。形の壊れた `practice` は初期値にする
  - 足りない入れ物（`bookmarks`・`predictions`・`misc`・`simPrefs`）を補う
  - 解答の記録・復習状態・解いた問題・メモ・予測・誤概念・設定・連続日数は変えない（`tests/store.test.mjs`・`tests/e2e.mjs` で、実際の形の v0.5 の記録 `tests/fixtures/state-v0.5.json` を使って確認）
- 各手順はコピーに対して動く純関数で、何度実行しても同じ。移行の途中で失敗したら本体には書かない。
- **v0.5 のアプリが v0.6 の記録を開いたとき**：v0.5 の `load()` は schema を変えずに使うので、v0.6 の項目（ブックマークなど）も残したまま保存する。v0.5 で「読み込む」をすると schema 4 になるが、v0.6 はそれを v0.5 の記録として移行し直す（冪等なので結果は同じ）。
  同期では、v0.5 のアプリは手元の schema（1）でサーバーに書き戻そうとするので、`guard.sql` を入れていれば CE003 で断られる（v0.6 の記録を v0.5 の形・5000件の切り捨てで上書きしない）。
- **決まり**：schema を 5 より上げる（v0.5 とも v0.6 とも違う形にする）ときは、保存先のキーも変える（`tests/store_edge.test.mjs` で検査）。
  手順：`SCHEMA` を上げる → `MIGRATIONS` に `{ from: 旧, up }` を足す → `KEY` を変え、起動時に旧キーから読む処理を足す → `tests/store.test.mjs` に旧版の記録を足す。

新しい schema の記録は、v0.6 以降のアプリは上書きしない（起動時に読んだ・保存しようとした・別のタブが書いた・サーバーにある、のどれでも）。

## 5. バックアップと復元

書き出し（設定の「書き出す」）の形（v1）：

```json
{ "format": "celestia-backup", "version": 1, "app": "0.6.0", "schema": 5,
  "createdAt": "2026-11-24T…Z", "counts": { "log": 60, "pstat": 24, "problems": 24, "bookmarks": 0, "explains": 2, "predictions": 2, "misc": 2 },
  "hash": "fnv1a32:xxxxxxxx", "state": { …celestia.v1 と同じ… } }
```

検査（`readBackup`。断るときは理由と `code`）は Lemma と同じ：`format`（`not_celestia`。Lemma のバックアップ `lemma-backup` も断る）・`version`・`state`・大きさ（`too_large`）・`app`・`schema`・`createdAt`・`hash`（`no_hash`/`bad_hash`/`hash_mismatch`）・`counts`（7項目。実際の件数との食い違いは注意）・値の点検（`invalid_values`）。
`hash` は FNV-1a（32bit）で、壊れた・途中で切れたファイルを見分けるためのもの（改ざん防止ではない）。

- 旧形式（v0.5 までの書き出し＝記録そのもの。schema 1 か 4）も読める（移行してから使う。hash による確認はできない）。
- v0.5 のアプリは v1 形式のファイルを読み込まない（`problems` と `log` が無いので「読み込めないファイル」になる。記録を壊さない）。
- 読み込むときは `settingsAt`・`resetAt`・`localResetAt` を捨て、設定（`settings`・`practice`・`simPrefs`）は今のものを残す。
- 「合わせる」か「置き換える」を選ぶ。置き換えでは `resetAt` を付ける。今の記録は控えに残し（残せなければファイルに書き出す）、保存に失敗したら元に戻す。件数の変化と点検の結果を表示する。
- 「記録をすべて消す」も控えを残す。「ログアウトしてこの端末だけ消す」「ログアウト時に端末の記録も消す」「アカウントの記録だけを使う」は、共用端末のため控えを残さず、同期できたことを確かめてから消し、`localResetAt` を付ける。

## 6. 同期（Supabase）

1回の同期：サーバーの行を取る → 確かめる → 手元と合わせる → 違えば `rev` つきで書き戻す（版がずれていれば取り直して最大4回）。
ログイン直後は、端末とアカウントの両方に記録があるときだけ「合わせる」か「アカウントの記録だけを使う（端末の分は書き出してから消す）」かを聞く。

### 6.1 クライアント（アプリ）が保証すること

Lemma と同じ（`tests/sync.test.mjs`・`tests/server_guard.test.mjs` の後半）：壊れた・新しい schema のサーバーの記録とは合わせず書き戻さない／手元の記録が壊れていれば送らない／条件つき PATCH／同時に呼んでも通信は1回分／通信中のログアウトで手元を変えない／`localResetAt` は送らない。

失敗の分類（`classify()`）：

| kind | 例 | 動き |
|---|---|---|
| offline | 通信できない | 15・30・60・120秒あけて自動でやり直す |
| conflict | 同時の書き込み・CE002 | 同上 |
| server | 5xx・429・408 | 同上 |
| auth | ログイン切れ・通信中のログアウト | やり直さない |
| rejected | サーバー（guard.sql）が断った CE001・CE003〜CE010 | やり直さず、自動の同期を止める。理由と対処を表示し、「今すぐ同期」で再開 |
| data | サーバーか手元の記録が壊れている・新しい版 | 同上 |

合わせ方（`mergeState`。どちらから合わせても同じ・何度合わせても同じ）：

| 項目 | 規則 |
|---|---|
| 解答の記録 | 和集合。`(p, t)` が同じものは1つにし、項目は足す。**切り捨てない**（v0.5 は5000件で切り捨てていた） |
| 復習状態 | 2台で別々に復習していれば、解答の記録から古い方の状態を `review()` し直す。そうでなければ新しい方（v0.5 は新しい方だけ） |
| 解いた問題 | 最後の時刻は新しい方、回数は多い方 |
| ブックマーク | 付け外しの新しい方（外した印も残す） |
| メモ | 和集合（問題ごとに新しい5件） |
| 予測の記録 | 和集合（`(p, t)`）。**切り捨てない**（v0.5 は2000件で切り捨てていた） |
| 誤概念 | 回数は多い方・最後の時刻は新しい方・ほかの項目（`cleared` など）は新しい方 |
| 設定・選択・シミュレーションの条件 | `settingsAt` が新しい方 |
| 連続日数 | 新しい日付の方 |
| `resetAt` | 新しい方。印の古い側からは、その時刻より前の記録（解答・復習状態・解いた問題・ブックマーク・メモ・予測・誤概念）を持ち込まない |
| `localResetAt` | 手元の値だけ |

### 6.2 サーバー（`supabase/guard.sql`）が保証すること

`schema.sql` だけでは、行単位のアクセス制限（本人の行だけ）と大きさ（5,000,000 バイト未満）しか無い。
`guard.sql` は表 `celestia_state` への書き込みごとに次を確かめ、破ればエラーで断る（表は変えない）。Lemma v0.26 の規則と同じ考え方：

| code | 断る書き込み |
|---|---|
| CE001 | `state` がオブジェクトでない・`log` が配列でない・解答の記録に文字列の `p` と数値の `t` が無い・pstat などがオブジェクトでない・`predictions` が配列でない／`p`・`t` の無い予測 |
| CE002 | 作成で `rev` が 1 でない・更新で `rev` が「前の rev + 1」でない |
| CE003 | `schema` が 1〜100 の整数でない・更新で下げる（v0.5 のアプリの書き戻しもこれで断る） |
| CE004 | 前にあった解答の記録（`(p, t)`）が無くなる |
| CE005 | 残っている解答の記録から項目を消す（値の更新は許す。`resetAt` より前の記録は対象外） |
| CE006 | `resetAt` が数値でない・戻す・進めるときにサーバーの時計で1日より先の未来 |
| CE007 | pstat・problems・misc（`last`）、bookmarks（`t`）の項目が消える・時刻が戻る、problems・misc の回数 `n` が減る、メモが消える（問題ごとに新しい5件を超えて押し出された古いメモは除く） |
| CE008 | 同じ `(p, t)` の解答の記録を増やす |
| CE009 | 同じ時刻のままの書き換え：pstat（値が範囲外・同じ `last` で reps を減らす・同じ `last`・`reps` で lapses を2以上変える）、bookmarks（同じ `t` で値を変える）、メモ（同じ `t` で本文を変える） |
| CE010 | 前にあった予測の記録（`(p, t)`）が無くなる |
| — | `localResetAt` は送られても保存しない |

`resetAt` を進めた書き込み（すべて消す・置き換え）だけ、その時刻より前のものは消えてよい。
守らないもの：settings・practice・simPrefs（最後に変えた方）、streak、古いメモ（新しい5件から押し出される）、誤概念の `cleared` などの状態。
守れないもの（Lemma と同じ）：`last` を進めたうえでの復習状態の値、新しく足された項目の中身、解答の記録の値の書き換え（再採点で正しく変わるため）。

**Lemma と共通のプロジェクトであることへの配慮**：`guard.sql` が作るのは schema `celestia_private`・その中の関数（`guard_version` と `celestia_` で始まる関数）・表 `celestia_state` のトリガー `celestia_state_guard` だけ。
Lemma の表 `lemma_state`・schema `lemma_private`・トリガー `lemma_state_guard`・関数 `public.delete_my_account()` には触れない（`tests/server_guard.test.mjs` で、Lemma の guard.sql を入れたデータベースに CELESTIA の SQL を2回入れ、Lemma の関数・トリガー・権限・制約・記録が1つも変わらず、Lemma の規則が今まで通り効くことを確認）。
エラーの code も別（Lemma は LM00x、CELESTIA は CE0xx）。

### 6.3 Supabase への適用手順

**実行する SQL と順番**（Supabase の SQL Editor。プロジェクトは Lemma と共通の「lemma」）：

1. `supabase/schema.sql` … 表 `celestia_state` を作る（すでにあれば作らない）。v0.5 の schema.sql で作った表は、制約 `state_size` を `celestia_state_size` に、ポリシーを `celestia_state: …` という名前に付け替える（中身は同じ。`celestia_state` の分だけ）。何度実行してもよい
2. `supabase/guard.sql` … サーバー側の書き込み制限。何度実行してもよい（関数は置き換え、トリガーは作り直す）。表と既存の記録は変えない

**適用済みかの確認**：

```sql
select celestia_private.guard_version();   -- '0.6.0' なら入っている（無ければ「関数が無い」エラー）
select t.tgname, s.nspname, p.proname from pg_trigger t join pg_proc p on p.oid = t.tgfoid join pg_namespace s on s.oid = p.pronamespace
where t.tgname in ('celestia_state_guard', 'lemma_state_guard');   -- celestia_state_guard | celestia_private と lemma_state_guard | lemma_private の2行
```

1. **適用前の確認**（何も変えない）：

   ```sql
   select user_id, rev, jsonb_typeof(state) as state_type, jsonb_typeof(state -> 'log') as log_type, state -> 'schema' as schema,
     case when jsonb_typeof(state -> 'log') = 'array' then jsonb_array_length(state -> 'log') end as log_n,
     (select count(*) from jsonb_array_elements(case when jsonb_typeof(state -> 'log') = 'array' then state -> 'log' else '[]'::jsonb end) e
        where jsonb_typeof(e) <> 'object' or jsonb_typeof(e -> 'p') is distinct from 'string' or jsonb_typeof(e -> 't') is distinct from 'number') as bad_entries,
     jsonb_typeof(state -> 'predictions') as pred_type, pg_column_size(state) as bytes
   from public.celestia_state order by bytes desc;
   ```

   `state_type` が `object`、`log_type` が `array`、`schema` が整数、`bad_entries` が 0、`pred_type` が `array` か空なら問題ない。

2. **念のための写し**（任意・推奨）：

   ```sql
   create table public.celestia_state_backup_v06 as select * from public.celestia_state;
   alter table public.celestia_state_backup_v06 enable row level security;   -- API から見えないように
   ```

3. **適用**：`supabase/schema.sql`、続けて `supabase/guard.sql` の中身をそれぞれ貼り付けて実行する。

4. **適用後の確認**：

   ```sql
   select celestia_private.guard_version();                                    -- '0.6.0'
   select lemma_private.guard_version();                                       -- Lemma の版（'0.26.0' など）が変わっていない
   select p.proname, has_function_privilege('anon', p.oid, 'EXECUTE') as anon, has_function_privilege('authenticated', p.oid, 'EXECUTE') as authenticated,
          p.proconfig, p.prosecdef from pg_proc p join pg_namespace s on s.oid = p.pronamespace where s.nspname = 'celestia_private';
   -- anon = false・authenticated = true・proconfig = {search_path=""}・prosecdef = false
   select conname from pg_constraint where conrelid = 'public.celestia_state'::regclass;   -- celestia_state_size・主キー・外部キー（celestia_state_user_id_fkey）
   select policyname from pg_policies where tablename = 'celestia_state';                    -- celestia_state: … の4つ
   select count(*) from public.celestia_state a full join public.celestia_state_backup_v06 b using (user_id)
   where a.state is distinct from b.state or a.rev is distinct from b.rev;                   -- 0（既存の記録が変わっていない）
   begin;   -- 規則が効いている（取り消す）
   update public.celestia_state set state = jsonb_set(state, '{log}', '[]'::jsonb), rev = rev + 1
   where user_id = (select user_id from public.celestia_state where jsonb_array_length(state -> 'log') > 0 limit 1);
   rollback;                                                   -- ERROR: celestia: 前にあった解答の記録 … （CE004）
   ```

   **PostgREST 経由で code が返るかの確認**（本番の自分の行で。記録はそのまま、版番号だけを2つ飛ばす要求）：アプリにログインし「今すぐ同期」を押してから、開発者ツールのコンソールで

   ```js
   const { SUPABASE_URL: U, SUPABASE_KEY: K } = await import('./js/config.js');
   const a = JSON.parse(localStorage.getItem('celestia.auth'));
   const h = { apikey: K, Authorization: `Bearer ${a.access_token}`, 'Content-Type': 'application/json' };
   const [row] = await (await fetch(`${U}/rest/v1/celestia_state?select=state,rev&user_id=eq.${a.user.id}`, { headers: h })).json();
   const r = await fetch(`${U}/rest/v1/celestia_state?user_id=eq.${a.user.id}&rev=eq.${row.rev}`,
     { method: 'PATCH', headers: { ...h, Prefer: 'return=representation' }, body: JSON.stringify({ state: row.state, rev: row.rev + 2 }) });
   console.log(r.status, await r.text());
   ```

   期待：`400 {"code":"CE002", …}`。`200` なら guard が効いていない（rev が2進むだけで記録は変わらない）。

5. **形の崩れた行がある場合**：利用者に書き出してもらい、写しを取ってから直す。トリガーを一時的に外すときは必ず戻す：
   `alter table public.celestia_state disable trigger celestia_state_guard;` … `alter table public.celestia_state enable trigger celestia_state_guard;`

**戻し方**：`drop trigger if exists celestia_state_guard on public.celestia_state;`（Lemma の制限には影響しない）

**適用後に同期できなくなる利用者**：v0.5 のアプリのままの端末（CE003・CE004。アプリを更新する＝再読み込み）、端末の時計が1日以上進んだ状態で「すべて消す」をした人（CE006。日時を直してもう一度）。

### 6.4 PGlite で確かめたこと・実際の Supabase でしか確かめられないこと

`tests/server_guard.test.mjs`（PostgreSQL の WASM 版 PGlite）：上の規則と通るべき書き込みのすべて、本物の `sync.js`・`store.js` によるトリガーつきの表での同期（2台・同時更新・2台の復習・障害と再試行・全消去・古い端末からの復活防止・置き換え・ログアウトして消した端末・5000件超・断られた後の停止と再開・v0.5 の書き戻し）、権限（anon は使えない・search_path が空・SECURITY DEFINER でない）、Lemma と共通のデータベースで Lemma 側が変わらないこと、`schema.sql` の名前の付け替えが何度実行しても同じで Lemma のポリシー・制約に触れないこと。
実際の Supabase でしか確かめられないこと：SQL を実行したこと、PostgREST がトリガーのエラーを 400 と `{ code: "CE0xx" }` で返すこと、RLS・認証と組み合わせた動き、本番の記録の大きさでの時間。

## 7. 同じ端末の複数のタブ

Lemma と同じ。保存のたびに、保存されている記録がこのタブの最後に読み書きしたものと違えば合わせてから書く（`mergeTabs`：`mergeState` に加えて `localResetAt` より前の記録を印の古い側から持ち込まない）。別のタブが保存すると `storage` イベントで合わせる。

## 8. 整合性の点検

`audit(state, { problemIds, now })` → `{ errors, warnings, stats }`。

| 区分 | code |
|---|---|
| エラー | `bad_container` `bad_log` `bad_predictions` `bad_schema` `log_bad_time` `log_bad_grade` `pstat_invalid` `problems_invalid` `bookmark_invalid` `explain_invalid` `misc_invalid` `reset_invalid` `settingsat_invalid` |
| 注意（auto） | `log_duplicate` `log_unsorted` `problems_behind_log` `prediction_duplicate` |
| 注意（human） | `log_unknown_problem` `pstat_unknown_problem` `bookmark_unknown_problem` `prediction_unknown_problem` `log_future` `reset_future` `streak_invalid` `schema_newer` |

`repair()` は削除しない（重複の統合・並べ替え・回数の引き上げだけ）。アプリは自動では直さない（設定に結果を表示するだけ）。

```
node tools/audit_state.mjs celestia-2026-11-24.json             # 点検（エラーがあれば終了コード 1）
node tools/audit_state.mjs celestia-2026-11-24.json --json
node tools/audit_state.mjs celestia-2026-11-24.json --fix 直した.json
```

## 9. 保存容量

Lemma と同じ判断：今は IndexedDB に移さない。本体が300万文字を超えたら1回知らせ、設定の「記録の点検」に保存している量を表示し、容量不足では控えを先に捨てる。
解答1件 約110文字＋問題ごとの復習状態・解いた回数 約120文字（947問で頭打ち、約11万文字）なので、1日10問でも10年以上余裕がある（localStorage は多くのブラウザで約500万文字。実機では確かめていない）。

## 10. 分かっている制約・残るリスク

- **v0.5 のアプリ**：同期で手元の schema（1）を書き戻し、解答の記録を5000件で切り捨てる。`guard.sql` を入れていれば断られる（CE003・CE004）が、入れていなければ通る。全端末を更新する必要がある。
- `guard.sql` は Supabase に入れるまで効かない。実際の Supabase（PostgREST 経由）では確かめていない。
- `resetAt` は端末の時計で付ける（時計が遅れた端末で消すと、ほかの端末の一部の記録が残る。進んだ時計は1日まではそのまま通る）。
- localStorage には複数タブにまたがるロックが無い（読み返しと次の保存で合わせ直すが、その前にタブを閉じるとその分は失われうる）。
- バックアップの `hash` は改ざん防止ではない。Safari は7日以上開かないとサイトの保存領域を消すことがある（ホーム画面に追加すると対象外）。
- 誤概念（`misc`）と予測（`predictions`）を書き込む画面はまだ無い（データの形と保護だけ先にある）。
