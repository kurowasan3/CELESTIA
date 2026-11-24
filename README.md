# CELESTIA 高校物理

<img src="icons/icon-192.png" width="96" alt="CELESTIA のアイコン：楕円軌道と、その焦点にある恒星">

図・シミュレーション・問題演習で、現象から理解する高校物理アプリ（試作版 v0.6）。
大学受験数学アプリ [Lemma](https://github.com/kurowasan3/lemma) の仕組み（復習・自動採点・同期・オフライン）を引き継いでいます（v0.6 で Lemma v0.28 の学習記録の保全・今日の学習・解答の形式などに揃えた。解法カードは引き継がない）。
iPhone・iPad・PC のブラウザで動き、ホーム画面に追加するとアプリとして使えます（PWA）。

目指す学び方：問題を読む → 図で現象をつかむ → シミュレーションで条件を変える → 結果を予測する → 法則で考える → 解く → 解答図と解説で深める → 間違えた概念を復習する

## いまできること

**v0.6 で Lemma v0.28 から引き継いだもの**

- **今日の学習**：ホームの一番上に、今日解いた数／目標、中断した演習の「前回の続き」、復習待ちを先に新しい問題で埋める「今日の演習を始める」（出題の理由つき。データが十分な苦手な単元から最大2問）。解答後は「類題を1問追加」で同じ単元・近い難易度・共通の概念や法則の問題を次に差し込める。自信なしで正解した問題は早めに復習に出す
- **解答の形式**：設定 →「解答の形式」で **記述式**（入力して自動採点）・**選択式**（選択肢から選ぶ。やり直しなし）・**自己採点**（紙に解いて答え合わせ）。選択肢のない欄は入力、記述・図示は自己採点。記録に形式も残る
- **判定がおかしい？**・**自信度と結果のずれへの一言**・**頻出度**（高・中・低。演習では未着手の問題のうち頻出度の高いものを優先）
- **ブックマーク**（問題画面の ☆、ホームから出題、同期で揃う）・**検索の並べ方と絞り込み**（関連度・頻出度・やさしい・難しい／すべて・未着手・間違えた・復習待ち・ブックマーク）
- **定着度**：得意・苦手（5問から判定、10問までは暫定）、最近の伸び、次にやること
- **学習記録の保全**：読めない記録は上書きせず退避、保存の失敗を知らせる、別のタブと合わせてから保存、移行・読み込み・消去の前に自動の控え、書き出し（形式・版・確認用の値つき）と読み込み（形を確かめてから「合わせる」か「置き換える」）、記録の点検、「記録をすべて消す」（ログイン中はアカウントも／この端末だけ）、ログイン時の確認（合わせる／アカウントの記録だけ）。v0.5 の記録はそのまま移行する。詳しくは [docs/data-integrity.md](docs/data-integrity.md)
- **サーバー側の書き込み制限**（`supabase/guard.sql`）：古い版のアプリなどが記録を消す・巻き戻す書き込みを、Supabase が断る

**v0.1〜v0.5**


- **図つきの問題**：問題文の途中に図（SVG）を置ける。1問に複数の図
- **解答図**：解答・解説にも図。答え合わせを押すまで表示せず、ファイルも読み込まない
- **答えを入力して自動採点**：Lemma の採点をそのまま使い、物理向けに `v_0`・`v0`（添字）と単位 `[m/s]` を追加
- **法則・誤概念**：問題は単元だけでなく「確かめる概念」「使う法則」「よくある誤解」につながる。「概念」タブは法則の一覧（分野別）で、法則のページから関連問題へ
- **全分野からランダム／単元ごと／復習／検索／問題セット**（Lemma と同じ）
- **自信度、定着度、ライト・ダーク、アカウントと同期、オフライン**（Lemma と同じ）
- **シミュレーション（全分野 24本）**：パラメータを動かすと、アニメーションとグラフが同じ物理モデルから一緒に変わる。再生・一時停止・巻き戻し（シークバーやグラフをなぞる）・速度変更、測定値つき
  - 力学（7）：斜方投射・斜面と摩擦・ばね振り子・単振り子・1次元の衝突・等速円運動・万有引力と軌道
  - 波動（5）：正弦波の伝わり方・弦の固有振動（定常波）・うなり・ドップラー効果・水面波の干渉
  - 光（3）：屈折と全反射・レンズの像・ヤングの実験
  - 電磁気（5）：点電荷の電場と電位・コンデンサーの充電と放電・磁場中の荷電粒子・磁場中を動く導体棒・交流回路（RLC 直列）
  - 熱（2）：気体分子の運動と圧力（壁の力積から圧力を「測る」）・気体の状態変化と熱力学第1法則
  - 原子（2）：光電効果・放射性崩壊と半減期
- **問題とシミュレーションの連携**：解いたあとに「シミュレーションで確かめる」で、問題と同じ条件のシミュレーションが開く。グラフ読解の問題のグラフ（v-t、波形 y-x など）も、シミュレーションのモデルから描く
- **法則の導出（68件）**：導出できる法則は1段ずつ開ける導出、出発点となる法則・実験から得た法則は根拠と由来。法則の一覧・解説から法則のページへ
- **問題 947問**：すべての単元（24単元）に問題がある。基礎〜入試標準は 444問（各単元に15問以上）。うち 568問（ID が PHY- で始まる、入試標準〜物理探究）は問題管理表から取り込んだもの。分野は力学・力学総合・熱力学・波動・電磁気・電磁気学総合・原子。誤概念を扱う概念問題、グラフ読解、融合問題を含む。誤概念 23件

難易度は 0=基礎・1=標準・2=入試基礎・3=入試標準・4=入試発展・5=物理探究 の6段階（Lemma と同じ考え方。共通テスト型かどうかは難易度ではなく `examStyle` で持つ）。
難易度とは別に、問題ごとに**頻出度**（高・中・低）がある。入試で繰り返し出る考え方・手法を扱う問題ほど高い（全947問を判定済み：高 172・中 267・低 508。基準は [docs/data-spec.md](docs/data-spec.md#頻出度)）。
入力して採点できる問題（368問）には**選択式の選択肢**（各欄に誤答3つ）がある（`content/choices.json`、1585欄）。

Lemma から外したもの：計算トレーニング、グラフ電卓（2D・3D）、数学の内蔵問題。計算問題は通常の問題形式（`format: "calc"`）として扱います。

## フォルダ

```
content/                 正規データ（マスター。ここを編集する）
  taxonomy.json          分野・単元（ID 固定）
  concepts.json  laws.json  misconceptions.json  simulations.json  tags.json  collections/
  choices.json           選択式の選択肢（tools/make_choices.mjs が作る）
  values.json            頻出度の判定（tools/values.py で点検・更新）
  problems/<分野>/<問題ID>/
    problem.json         問題・答え・解説・図の情報
    q-*.fig.json         作図仕様（tools/figgen.py が SVG にする）
    q-*.svg  a-*.svg     図（q- は問題用、a- は解答用。名前は自由だが分けておくと見やすい）
data/                    アプリ用（生成物。tools/build_data.py が作る。手で編集しない）
js/                      アプリ（app.js・catalog.js・answer.js・parser.js・problem/figures.js・sync.js …）
  sim/core/              シミュレーションの共通部分：時間・パラメータ（単位変換）・軌跡・数値積分・イベント検出
  sim/models/mechanics/  物理モデル（画面に依存しない。Node だけでテストできる）
  sim/render/  sim/plot/  sim/ui/   描画・連動グラフ・プレーヤー
tools/                   変換・作図・検査（Python と Node）
tests/                   検査（node tests/run.mjs。--e2e でブラウザの通し操作も。npm install --no-save --no-package-lock @electric-sql/pglite playwright が必要）
docs/                    データ仕様（data-spec.md）・設計（DESIGN.md）
supabase/schema.sql      同期用の表（Supabase の SQL Editor で実行）
supabase/guard.sql       サーバー側の書き込み制限（schema.sql のあとに実行）
privacy.html terms.html  プライバシーポリシー・利用規約
js/store.js              学習記録の保存・移行・点検・バックアップ・合わせ方
```

## 問題を追加・修正する

```sh
pip3 install jsonschema
python3 tools/figgen.py                 # 作図仕様（*.fig.json）から SVG を作る
node tools/make_choices.mjs             # 入力用の答えから選択式の選択肢を作り直す（content/choices.json）
python3 tools/values.py                 # 頻出度の判定が全問そろっているか（無い・古いものを知らせる）
python3 tools/build_data.py             # content/ → data/
node tests/run.mjs                      # 検査（全部 ✓ になればOK）
python3 tools/validate_content.py --write-golden   # 新しい問題ID を台帳に足す（問題を足したとき）
```

1. `content/problems/<分野>/<新しい問題ID>/problem.json` を、既存の問題をまねて作る（形は [docs/data-spec.md](docs/data-spec.md)）
2. 図は作図仕様（`*.fig.json`）で書くのがおすすめ。力の向き・斜面の角度などを数値で書くので、図と物理条件の照合が自動になる
3. 上の行を実行。`validate_content.py` が、図の参照切れ・解答図の漏れ・図と物理条件の食い違い（力のつり合い・ベクトルの分解・角度）などを見つける
4. `content`・`data`・`tests` フォルダごとリポジトリにアップロード

規則：問題ID は変えない・消さない（消すと台帳 `tests/golden/ids.txt` との照合で検査が失敗する）。出さない問題は `quality.status` を `hidden` に。

### 選択式の選択肢（content/choices.json）

選択肢は「入力用の答え」（`answer.parts[].grading`）から `tools/make_choices.mjs` が決まった手順で作ります（AI ではないので、何度実行しても同じ）。Lemma の同名のツールの物理版で、候補は次の順に使います。

1. 問題に書いてある典型的な誤答（`answer.parts[].wrong`。誤概念に結びついたもの）
2. 式の形の誤り：sin と cos の取り違え・分母と分子の逆転・2倍や半分の付け忘れ（項ごとにも）・2乗の付け忘れ・和と差の取り違え・符号（仕事・電位・変化量など符号に意味がある量、負の答えだけ）
3. 数の誤り：単位つきの数値は 2倍・半分・指数の1つずれ（10倍）などを、有効数字の桁数をそろえて。回・個・倍などの小さい整数は ±1 も。分数は分母・分子の2倍や逆数。**単位はそのまま残す**

アプリの採点（`js/answer.js` の `checkField`）で不正解になること・候補どうしが表示も値も重ならないこと・数値は互いに 10% 以上離れていること・√4 や 2/4 や 6E/4 のような不自然な形でないことを確かめて、各欄に誤答を3つ選びます（3つ目は誤答をさらに1か所変えたものにして、正解だけが「真ん中」に見えないようにする）。2つしか作れない欄は選択肢3つ、1つ以下の欄は選択肢なし（その欄は入力式）。

```sh
node tools/make_choices.mjs             # content/choices.json を作り直す
node tools/check_choices.mjs            # 誤答がすべて不正解・正解が正解になることを採点で確かめる
python3 tools/build_data.py             # data/ に反映（本文の ch、一覧の ch: 1）
```

入力用の答えを直した問題は、`choices.json` の `spec` が合わなくなるので選択肢が出なくなります（`check_choices.mjs` が知らせる）。上の3行で作り直してください。

### 頻出度（content/values.json）

問題ごとの頻出度（1=低・2=中・3=高）は、問題の単元・テーマ・タグ・法則・問題文の中心の考え方を見て、大学入試でどれだけ繰り返し出るかを判定したものです（基準は [docs/data-spec.md](docs/data-spec.md#頻出度)）。判定には問題文のハッシュが付いていて、問題文を直すとその判定は使われなくなります（中の扱い）。

```sh
python3 tools/values.py                 # 点検（判定が無い・古い問題があれば知らせる。検査でも使う）
python3 tools/values.py --stale         # 判定し直す問題の一覧（テーマ・タグ・法則つき）
python3 tools/values.py --set PHY-IND-001=3 MEC-ENE-105=1   # 判定を書く（今の問題文のハッシュを付ける）
python3 tools/values.py --stats         # 分野ごとの分布
```

問題ごとに決めたいときは `problem.json` に `"value": 1〜3` を書けば、それが判定より優先されます。

## シミュレーションを追加する

1. `js/sim/models/<分野>/<名前>.js` に物理モデルを書く（SI 単位。解析解があれば `at(p, t)`、なければ `init`・`accel`（保存系は速度ベルレ）または `deriv`（RK4）。着地などは `event`）
2. `js/sim/render/` に描画器（`bounds` と `draw`）を書き、`RENDERERS` に登録する
3. `content/simulations.json` にパラメータ（単位・範囲・既定値）・量・グラフ・表示名を登録し、`status` を `ready` に
4. `tests/sim.test.mjs` に解析解・保存則との比較を足す（登録簿との整合は自動で検査される）
5. `sw.js` の `FILES` にファイルを足す

モデルの数値計算の精度は `tests/sim.test.mjs` で確かめています（例：振れ角 90° の単振り子の周期の比 1.18034、50周期でのエネルギーのずれ 1e-5 未満、楕円軌道での角運動量の保存、ドップラー効果で波面が届く間隔から測った振動数＝式、気体分子の壁への力積から測った圧力≒Nm⟨v²⟩/3V、RC 回路の熱＝½CE²、導体棒の F の仕事＝運動エネルギー＋熱、P-V 図の面積＝仕事）。

時間によらないシミュレーション（レンズ・屈折・干渉縞・電場・光電効果）は、モデルに `static: true` と `profile(p)` を書く。波形のように時刻ごとに形が変わるグラフは `snapshot(p, t)` を書き、登録簿のグラフに `"snap": true` を付ける。

## アカウントと同期（Supabase）

- Supabase プロジェクトは **Lemma と共通**（`js/config.js` の値も同じ）。アカウントも共通で、記録は別の表 `celestia_state` に保存
- 準備：`supabase/schema.sql`、続けて `supabase/guard.sql` を SQL Editor で実行（どちらも何度実行してもよい。名前はすべて `celestia_` で始まり、Lemma の表・関数・トリガーには触れない）。確認と戻し方は [docs/data-integrity.md](docs/data-integrity.md) §6.3
- 退会はアカウントごと消すので、Lemma の記録も消える（画面でもそう伝えている）
- 端末の保存キーは `celestia.v1`（記録）・`celestia.snap.v1`（控え）・`celestia.session.v1`（中断した演習）・`celestia.auth`。同じドメインの Lemma（`lemma.v2`）・ぶつりば（`lemma-physics.v1`）と混ざらない
- サーバーに断られたときのエラーは CE001〜CE010（Lemma の LM00x とは別）

## 配信用の dist/（Cloudflare Pages）

Lemma と同じく、アプリ本体と data/ だけを集めた `dist/` を配信できます（リポジトリを非公開にしても公開できる）。

```sh
python3 tools/build_dist.py             # dist/ を作り直す（出力先を変えるときは python3 tools/build_dist.py OUT）
node tests/e2e_dist.mjs                 # dist/ を CSP つきで配信して、採点・KaTeX・図・シミュレーションが動くか確かめる
```

- 入れるもの：`index.html`・`privacy.html`・`terms.html`（あれば）・`manifest.webmanifest`・`sw.js`・`css/`・`js/`（`sim/`・`problem/` を含む）・`vendor/`・`icons/`（作業用の `icons/src` は除く）・`data/*.json`・`data/assets/`（問題の図）
- 足すもの：`_headers`（CSP などのセキュリティヘッダー。Cloudflare Pages・Netlify の書式）・`robots.txt`
- 入れないもの：`content/`・`tools/`・`tests/`・`docs/`・`supabase/`・`verify/`・`README.md`（入っていたら build_dist.py がエラーにする）
- `.github/workflows/test.yml` は、main に push して検査が通ると `dist/` を Cloudflare Pages に配置する。リポジトリの Settings → Secrets and variables → Actions で Variables `CF_PAGES_PROJECT`、Secrets `CLOUDFLARE_API_TOKEN`・`CLOUDFLARE_ACCOUNT_ID` を登録するまでは、検査だけが動く（検査の結果の dist は Actions の成果物としても7日間取り出せる）

## GitHub Pages で公開する

Lemma と同じ手順です。リポジトリ（例：`mythos0217/celestia`、Public）の直下にこのフォルダの中身を置き、Settings → Pages で `main` / `/ (root)` を選ぶと `https://mythos0217.github.io/celestia/` で開けます。

- `.github/workflows/test.yml`（自動検査）は、Web の「Upload files」でフォルダごとドラッグすると落ちることがあります。Actions タブに何も出ないときは、このファイルを個別にアップロードしてください
- 更新したときは `sw.js` の `VERSION` と、各 JS の `FILE_VERSION`、`css/style.css` の `--celestia-version` を同じ値に上げる
- オフライン用のキャッシュ名は `celestia-` で始まり、消すのも `celestia-` のものだけです（同じドメインの Lemma・ぶつりばのキャッシュを消さない）

## 問題管理表からの取り込み

問題管理表の形式（1問1フォルダの `item.json`・`problem.svg`・`solution.svg`）で作った問題は、次のコマンドで content/ に取り込める。

```
python3 tools/import_items.py <content フォルダ> --sheet CELESTIA問題管理表.xlsx --verify <verify フォルダ>
python3 tools/build_data.py && python3 tools/validate_content.py --write-golden && node tests/run.mjs --e2e
```

- 問題ID は item.json の id をそのまま使う（同じ ID は上書き）。単元は単元名で決まる。
- 難易度 1〜5 は level 基礎・標準・入試標準・入試発展・物理探究に対応させる（2 は標準、3 は入試標準）。
- 問題文の enumerate が小問、答えの「(1) … (2) …」が小問ごとの答えになる。入力採点はしない（自己採点）。
- 参考元はアプリに出さないので取り込まない。管理表の作成日・作成方法・検証3項目は provenance と quality に写す。
- 概念はタグ・テーマ・核の語から自動で付ける（最大3つ）。法則への結びつけはしない。
- 検証スクリプトは `verify/<ID>.py` に置く（アプリでは使わない）。
- 色を直接指定した SVG は、ダークモードでも白い紙の上に表示する。
