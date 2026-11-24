# CELESTIA データ仕様（v0.4）

正規データ `content/`（マスター）→ `tools/build_data.py` → アプリ用データ `data/`。アプリが読むのは `data/` だけ。
JSON Schema は `tools/schema/`、Schema で見られない関係の検査は `tools/validate_content.py`。

## ID の規則

| ID | 形 | 規則 |
| --- | --- | --- |
| 問題ID | `MEC-INC-001` | `<分野3字>-<単元略号3字>-<番号>`。一度使ったら変えない・再利用しない。分野を移しても ID はそのまま（フォルダだけ移す） |
| 分野ID | `MEC` `MSY` `THM` `WAV` `ELM` `ESY` `ATM` | 力学・力学総合・熱力学・波動・電磁気・電磁気学総合・原子（`SYN`「分野横断」は v0.4 で廃止し、問題は総合の2分野へ移した。ID は変えていない） |
| 単元ID | `MEC-02` | `taxonomy.json` で固定。使わなくなったら `retired` に残し、再利用しない |
| 概念ID・法則ID | `normal-force` `newton-2nd` | 英小文字とハイフン |
| 誤概念ID | `MC-NORMAL-EQ-MG` | `MC-` で始まる英大文字 |
| シミュレーションID | `mech.projectile` | `<分野>.<名前>` |
| 図ID | `q-setup` `a-forces` | 問題の中で一意。英小文字とハイフン |

## 問題（`content/problems/<分野>/<問題ID>/problem.json`、`celestia.problem/1`）

| 項目 | 必須 | 内容 |
| --- | --- | --- |
| `id` `schema` `field` `unit` | ○ | 問題ID、`"celestia.problem/1"`、分野ID、単元ID |
| `topic` `course` `minutes` | | テーマ、`物理基礎`／`物理`、想定時間（分） |
| `level` | ○ | 0=基礎 1=標準 2=入試基礎 3=入試標準 4=入試発展 5=物理探究 |
| `value` | | 頻出度 1=低 2=中 3=高（下の「頻出度」）。書けば `content/values.json` の判定より優先。ふつうは書かない |
| `format` | ○ | `choice` `numeric` `calc` `written` `derivation` `graph-read` `figure-read` `construct` `experiment` `concept` `synthesis`（複数可） |
| `examStyle` | | `common-test` `national-2nd` `private` `original` |
| `statement.blocks` | ○ | 本文（下の「ブロック」） |
| `figures` | ○ | 図の一覧（下の「図」）。無ければ `[]` |
| `answer.parts` | ○ | 小問ごとの答え `{id, tex, grading?, wrong?}`。小問が無い問題は `id: "main"` 1つ。`wrong`（`[{value, misconception?}]`、典型的な誤答を `grading` と同じ書式で）は選択式の選択肢の誤答に優先して使う |
| `solutions` | ○ | `{id: "main"|"alt1"|"note", kind, blocks}`。別解・補足を足せる |
| `idea` | | 核となる考え方 `{tex}` |
| `simulation` | | `{id, preset, locked, show}`。`show` は `before`／`after-attempt`／`solution` |
| `predict` | | 予測 → 実験 → 説明 `{prompt, choices[{id, tex, misconception?}], answer, experiment{sim, vary, measure}, explain{laws, tex}}` |
| `choices` | | 選択問題の選択肢。誤答の選択肢に `misconception` を付けると、選んだときに誤概念として記録できる（Phase 5） |
| `concepts` `laws` `misconceptions` | ○ | 概念・法則・誤概念の ID（無ければ `[]`） |
| `tags` `provenance` `quality` | ○ | `quality.status` は `draft`／`review`／`published`／`hidden`。アプリに出るのは `published` だけ |

### ブロック

```json
{"type": "text", "tex": "文章（LaTeX。$...$ が数式）"}
{"type": "figure", "ref": "図ID", "placement": "block | float-right | side-by-side"}
{"type": "parts", "items": [{"id": "1", "tex": "…", "figures": ["図ID"]}]}
{"type": "simulation", "caption": "…"}
```

並びがそのまま表示順。図の位置も、ブロックの位置で決まる（`order` を手で管理しない）。

### 入力用の答え（`answer.parts[].grading`）

Lemma の書式（`js/answer.js` の先頭）に、次の2つを足した（また、Lemma で積分定数として特別扱いしていた `C` は、物理では電気容量などのふつうの文字として扱う）：

- 添字：`v_0`、`v_{max}`、`m_A`。文字の直後の数字も添字とみなす（`v0` = `v_0`）
- 単位：末尾の `[ ]`。`6 [Ω]`、`v_0 cos(θ) [m/s]`（今は表示用。単位表は `tools/validate_content.py` の `UNITS`）

証明・説明など自動で採点しない小問は `grading` を書かない（自己採点になる）。

## 選択式の選択肢（`content/choices.json`）

`tools/make_choices.mjs` が入力用の答えから作る（手で編集しない。作り方は README）。Lemma と同じ形に、検査用の `t` を足したもの。

```json
{"ELM-CAP-001": {"spec": "(1) Q：ε_0 S V/d\n(2) (1/2)倍…",
                 "f": [["\\frac{\\varepsilon_{0}SV}{d}", "\\frac{2\\varepsilon_{0}SV}{d}", "…", "…"], …],
                 "t": [["ε_0 S V/d", "2ε_0 S V/d", "d/(ε_0 S V)", "6ε_0 S V/d"], …]}}
```

| 項目 | 内容 |
| --- | --- |
| `spec` | 作ったときの入力用の答え（`build_data.py` の `grading_spec()` と同じ文字列）。今の答えと違えば古い選択肢で、`data/` に出さない |
| `f` | 入力欄（全小問の欄を順に並べたもの）ごとに `[正解の TeX, 誤答の TeX, …]`（3〜4個）か `null`（選択肢なし＝その欄は入力式）。数値は書いたとおりの桁数、`×10^n`、単位は `\mathrm{}` |
| `t` | `f` と同じ並びの、入力用の答えの書式の文字列（単位つき）。`tools/check_choices.mjs` が、`t[0]` が正解・`t[1..]` がすべて不正解・互いに等しくないこと、`t` から作った TeX が `f` と同じこと、KaTeX で描けることを確かめる |

## 頻出度

`content/values.json`：`{問題ID: {"h": 問題文のハッシュ, "v": 1〜3}}`（Lemma と同じ形）。`h` は `statement` を JSON にしたものの SHA-1 の先頭12桁（`content_lib.statement_hash`）で、問題文を直すとその判定は使われなくなる（`tools/values.py` が「古い」と知らせる）。`problem.json` の `value` があればそれが優先。

判定は、問題の単元・テーマ（`topic`）・タグ・法則・問題文の中心の考え方を見て、**その考え方・手法が大学入試（共通テスト・国公立2次・私大）でどれだけ繰り返し出るか**で決める。難易度とは別の軸（基礎でも頻出度が低いもの、物理探究でも頻出度が高いものがある）。

| 頻出度 | 基準 | 例 |
| --- | --- | --- |
| 3 高 | 多くの大学で毎年のように出る定番の考え方・手法。問題の中心がそれそのもの | 曲面上の力学的エネルギー保存・鉛直面内の円運動・反発係数と運動量保存・台と小物体の運動量保存・ばね振り子・人工衛星とケプラーの法則・熱力学第1法則と P–V 図のサイクル・分子運動論の圧力の導出・ドップラー効果・ヤングの実験・薄膜干渉・気柱と弦の共鳴・コンデンサーのつなぎかえと誘電体・孤立部分の電荷保存・キルヒホッフの法則・導体棒の電磁誘導と終端速度・磁場中の荷電粒子の円運動・RLC と電気振動・光電効果・ボーア模型・質量欠損と半減期 |
| 2 中 | 入試に出るが3ほど多くない話題、または定番の考え方を見慣れない設定で使うもの | らせん運動・サイクロイド・ホール効果の応用・抵抗の測定と計器・電位差計・LC 回路のダイオード・自由電子モデルの拡張・二重スリットの物質波・地球を貫くトンネル・連成振動・断熱変化での気体の振動・2 物体の重心系での運動 |
| 1 低 | 特殊な装置・研究的な発展・一度きりの設定で、同じ考え方が入試でくり返し出るとは言いにくいもの。単発の知識（線量の単位など） | レーザー冷却・ゼーマン減速器・中性子干渉計・ポジトロニウムの3光子消滅・パラメトリック増幅・チャージポンプ・鏡像法・回転系の細かな安定性・表面張力の細部・多モード光導波路 |

アプリでは、問題画面に表示し、演習ではまだ解いていない問題のうち頻出度の高いものを優先して出す（Lemma と同じ）。判定は人の物理の判断によるもの（自動の判定ではない）。新しい問題を足したら `python3 tools/values.py --stale` で一覧し、判定して `--set` で書く。

## 図

```json
{"id": "a-forces", "role": "answer", "kind": "vector", "format": "svg", "src": "a-forces.svg",
 "alt": "物体にはたらく重力・垂直抗力・静止摩擦力", "caption": "…",
 "physics": {"angle_deg": 30, "equilibrium": true,
   "forces": [{"name": "mg", "dir_deg": 270}],
   "decompositions": [{"of": "force:mg", "into": ["component:mg-parallel", "component:mg-normal"]}]}}
```

| 項目 | 内容 |
| --- | --- |
| `role` | `question`（問題用）／`answer`（解答用）／`both`。**解答用の図は問題文から参照できない**（検査でエラー）。アプリは答え合わせまで表示もファイルの読み込みもしない |
| `kind` | `diagram` `vector` `circuit` `ray` `waveform` `apparatus` `state` `graph` `photo` `sim-initial` |
| `format` | `svg`（標準）、`png`／`webp`（写真）、`plot-json`（データから描くグラフ。未実装）、`sim-plot`（シミュレーションのモデルから描くグラフ：`"sim": {"id", "preset", "t"?}, "plot": {"x": "t", "y": ["vy"], "snap"?}`。波形 y-x は `"x": "x", "snap": true` と時刻 `t`）、`sim-snapshot`（未実装） |
| `src` | 問題フォルダ内のファイル名。Base64 の埋め込みは不可 |
| `alt` | 代替テキスト（必須。図が見えない人が問題を解けるように） |
| `physics` | この図が表す物理条件。`validate_content.py` が SVG の実際の形と照合する |

### 図と物理条件の照合（`physics`）

| 宣言 | 照合すること |
| --- | --- |
| `angle_deg` | `data-phys="surface"` の線の傾き（±2°） |
| `forces` `vectors` | `data-phys="force:N"` などの矢印の向き（±2°） |
| `equilibrium: true` | `force:` の矢印の和が 0（最大の力の 4% 以内） |
| `netForce: {dir_deg}` | 運動している物体：`force:` の矢印の和（合力）がその向き（±5°）で、0 ではない |
| `torque: {about: [x, y]}` | 剛体：その点のまわりの力のモーメント（矢印の始点＝作用点）の和が 0 |
| `refraction: {n1, n2, in, out}` | 光線：境界は水平として、入射光線と屈折光線の向きが屈折の法則 n₁sinθ₁ = n₂sinθ₂ を満たす（2% 以内） |
| `decompositions` | 成分の矢印の和が元の矢印に一致 |
| `equalLength` | 同じ大きさのはずの矢印の長さがそろう |
| `junctions` | 分岐点で電流の和が保存される |

見た目が正しくても物理的に誤った図（AI が作った図に多い）を止めるための仕組み。自動の照合に通っても `quality.checks.figures` は `"auto"` 止まりで、人が確認したら `"human"`。

### 作図仕様（`*.fig.json`、`celestia.figure/1`）

`tools/figgen.py` が SVG を作る。部品：`line` `poly` `arrow` `block` `ball` `arc` `text` `spring` `curve`（放物線）`resistor` `battery` `wire` `dot` `circle` `sine`（正弦波の曲線）。
座標は SVG と同じ（右が x、下が y）、角度は数学と同じ（右 0°、反時計回り正）。詳しくは `tools/figgen.py` の先頭。
`node tests/run.mjs` は、SVG が仕様から作り直したものと同じか（作り忘れ）も確かめる。

### 手描き・AI 生成の SVG の決まり

- `viewBox` を付ける。`script`・`foreignObject`・`on*` 属性・外部参照・埋め込み画像は使わない（検査でエラー）
- 色は `currentColor`（明暗の設定に合わせて変わる）。力・速度・成分・電流は `<g class="f-force|f-vel|f-comp|f-cur">`
- 照合したい矢印は `<g data-phys="force:N">` の中に、軸の `<line>` と、先端を最初の点にした `<polygon>` を置く
- 面は `<line data-phys="surface">`

## 概念・法則・誤概念・シミュレーション

| ファイル | 1件の形 |
| --- | --- |
| `concepts.json` | `{id, name, field, requires: [前提の概念], laws: [法則]}`（前提関係の循環は検査でエラー） |
| `laws.json` | `{id, name, tex, conditions: 適用条件, concepts, origin, basis?, derivation?, simulations?}`。`origin` は `fundamental`（出発点となる法則）／`empirical`（実験から得た法則）／`derived`（導出できる）／`definition`。`derived` は `derivation`（`[{tex}]`、1段ずつ）が必須、`fundamental`・`empirical` は `basis`（根拠・由来）が必須 |
| `misconceptions.json` | `{id, wrong: 誤った考え, right: 正しい考え, concepts, remedy: [補強問題の ID]}` |
| `simulations.json` | `{id, name, field, model: "mechanics/projectile", status: planned|ready, static?, params{名前: {label, unit, min, max, step, default}}, observables{量: 単位}, labels{量: 表示名}, measures{名前: {label, unit}}, plots}`。`plots` の各段は `[x, y]` か `{x, y | ys: [..], snap?, mark?, yRange?, label?}`（`ys` で重ね描き、`snap` で時刻ごとの分布、`mark` で目印の縦線）。`ready` なら `js/sim/models/<model>.js` が必要（検査） |

## アプリ用データ（`data/`、`celestia.data/1`）

| ファイル | 中身 |
| --- | --- |
| `index.json` | `format` `dataVersion` `levels`（`{"0": "基礎", …, "5": "物理探究"}`）`subjects`（分野の並び）`fields`（単元）`problems`（一覧：`id f lv t th tg fm c m sim pr in ch v ok`。`ch: 1` は選択式の選択肢がある問題、`v` は頻出度で 1 か 3 だけを書く（2＝中は省く）） |
| `f-<単元ID>.json` | 本文：`q a s k`（HTML）、`figs`（図の情報）、`in`（入力用の答え）、`ch`（選択式の選択肢。`choices.json` の `f` そのまま：入力欄の順に `[正解の TeX, 誤答の TeX, …]` か `null`。先頭が正解で、表示のときに混ぜる）、`sim` `predict` `choices`（選択問題の選択肢。`ch` とは別）`laws` `alt`（別解） |
| `knowledge.json` | 概念・法則・誤概念・シミュレーション登録簿 |
| `search.json` `collections.json` | 検索の索引（概念名・法則名も入る）、問題セット |
| `assets/<問題ID>/` | 図のファイル。名前に中身のハッシュ（`q-setup-35755270.svg`）が付くので、図を直すとキャッシュも入れ替わる |

`dataVersion` は content/ のファイル（`choices.json`・`values.json` を含む）の中身から計算する。選択肢や頻出度だけを直しても版が変わる。

本文の HTML では、図は `<figure class="pfig" data-fig="図ID"></figure>` の空の枠で、`js/problem/figures.js` が中身を入れる。

## 版は3つを分ける（Lemma と同じ）

| 版 | 場所 | いつ上げるか |
| --- | --- | --- |
| コードの版 | 各 JS の `FILE_VERSION`、CSS の `--celestia-version`、`sw.js` の `VERSION` | アプリを更新したとき（食い違うと赤い警告） |
| データ形式の版 | `index.json` の `format` | `data/` の項目の意味を変えたとき |
| 学習記録の版 | localStorage `celestia.v1` の `schema`（v0.6 から 5。v0.5 までは 1、v0.5 の読み込みでは 4） | 学習記録の形を変えたとき（`js/store.js` の `SCHEMA`・`MIGRATIONS`。5 より上げるときはキーも変える） |

## アプリが読む項目（v0.6）

- `index.json` の各問題の `v`（頻出度 1〜3。無ければ 2＝中）と `ch`（1 なら選択肢あり）。どちらも無くてもアプリは動く（`js/catalog.js` の `valueOf`・`getProblem`）
- 本文の `ch`：入力欄（`in` を `parseSpec` で読んだ欄の並び）と同じ長さの配列。各要素は `[正解, 誤答, …]`（TeX の文字列。先頭が正解、表示の順はアプリが混ぜる）か `null`（その欄は入力で答える）。Lemma v0.28 と同じ形
- 学習記録（localStorage `celestia.v1`・同期表 `celestia_state`）の形・移行・バックアップ・同期の規則は [data-integrity.md](data-integrity.md)
