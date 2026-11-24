// CELESTIA の検査をまとめて実行する（GitHub Actions でも同じものを実行）
//   node tests/run.mjs          データと単体テスト
//   node tests/run.mjs --e2e    ブラウザでの通し操作も（Playwright が必要）
// 必要なもの：npm install --no-save --no-package-lock @electric-sql/pglite playwright（server_guard と e2e）
// どれか1つでも失敗したら終了コード 1。
import { spawnSync } from 'child_process';
import fs from 'fs';
import os from 'os';
import path from 'path';
import { fileURLToPath } from 'url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const py = process.env.PYTHON || 'python3';
const tmp = fs.mkdtempSync(path.join(os.tmpdir(), 'celestia-'));

const steps = [
  ['正規データの Schema（問題・作図仕様・概念・法則・誤概念・シミュレーション）', py, ['tools/check_schema.py']],
  ['正規データの整合（ID の台帳・図の参照と漏れ・SVG の安全性・図と物理条件の照合・参照・単位）', py, ['tools/validate_content.py']],
  ['SVG が作図仕様（*.fig.json）から作り直したものと同じ', py, ['tools/figgen.py', '--check']],
  ['content/ から data/ を作り直すと同じになる', py, ['-c', `
import filecmp, sys
from pathlib import Path
sys.path.insert(0, 'tools')
import build_data as BD
t = Path(${JSON.stringify(path.join(tmp, 'data'))})
BD.build('content', t)
def files(d):
    return sorted(str(p.relative_to(d)) for p in Path(d).rglob('*') if p.is_file())
a, b = files('data'), files(t)
ok = a == b
if not ok:
    print('ファイルの違い：', sorted(set(a) ^ set(b))[:10])
for f in a:
    if f in b and not filecmp.cmp(Path('data') / f, t / f, shallow=False):
        print('中身が違う：', f); ok = False
if not ok:
    print('→ python3 tools/build_data.py を実行して data/ を作り直してください')
sys.exit(0 if ok else 1)`]],
  ['入力用の答え（読める・正解を入れると正解になる）', 'node', ['tools/check_answers.mjs', 'data']],
  ['すべての数式を KaTeX で描画できる', 'node', ['tools/check_katex.cjs']],
  ['自動採点の回帰テスト', 'node', ['tests/answer.test.mjs']],
  ['選択式の選択肢（content/choices.json）が入力用の答えから作り直したものと同じ', 'node', ['tools/make_choices.mjs', 'content', '--check']],
  ['選択式の選択肢を採点で確かめる（誤答はすべて不正解・正解は正解・重なりなし・KaTeX）', 'node', ['tools/check_choices.mjs']],
  ['頻出度（content/values.json）が全問・今の問題文に対してそろっている', py, ['tools/values.py']],
  ['配信用の dist/（必要なファイルだけ・CSP・配信してはいけないものが無い）', py, ['tools/build_dist.py', path.join(tmp, 'dist')]],
  ['catalog.js の単体テスト（共通ビュー・図・概念・検索・問題セット）', 'node', ['tests/catalog.test.mjs']],
  ['同期の単体テスト（記録の合わせ方・にせ Supabase との通信）', 'node', ['tests/sync.test.mjs']],
  ['すべての JS モジュールを読み込める（構文エラーがない）', 'node', ['tests/syntax.test.mjs']],
  ['シミュレーションの物理モデル（解析解・保存則・登録簿・依存の向き）', 'node', ['tests/sim.test.mjs']],
  // v0.6：Lemma v0.28 の学習記録の保全・推薦を引き継いだもの
  ['今日の演習の推薦（架空の学習者・境界条件）', 'node', ['tests/recommend.test.mjs']],
  ['学習記録の保存・点検・バックアップ・移行（store.js。v0.5 の記録の移行を含む）', 'node', ['tests/store.test.mjs']],
  ['学習記録の保全の境界条件（新しい schema・バックアップ形式・保存先の障害・複数タブ・移行の再試行）', 'node', ['tests/store_edge.test.mjs']],
  ['学習記録の点検ツール（v0.5 の実際の形の記録にエラーが無い）', 'node', ['tools/audit_state.mjs', 'tests/fixtures/state-v0.5.json']],
  ['サーバー側の書き込み制限（supabase/guard.sql を PostgreSQL＝PGlite で。Lemma と共通のプロジェクトでも Lemma 側を変えない）', 'node', ['tests/server_guard.test.mjs']],
];
if (process.argv.includes('--e2e')) steps.push(['ブラウザでの通し操作', 'node', ['tests/e2e.mjs']], ['解答の形式（記述式・選択式・自己採点）', 'node', ['tests/e2e_modes.mjs']], ['配信用の dist/ を CSP つきで動かす', 'node', ['tests/e2e_dist.mjs']]);

let failed = 0;
for (const [name, cmd, args] of steps) {
  const t0 = Date.now();
  const r = spawnSync(cmd, args, { cwd: root, encoding: 'utf8' });
  const ok = r.status === 0;
  const out = ((r.stdout || '') + (r.stderr || '')).split('\n')
    .filter(l => l.trim() && !/No character metrics|UserWarning|warn\(msg\)/.test(l));
  console.log(`${ok ? '✓' : '✗'} ${name}（${((Date.now() - t0) / 1000).toFixed(1)}秒）`);
  if (!ok) {
    failed++;
    console.log(out.slice(-25).map(l => '    ' + l).join('\n') || `    ${r.error || '終了コード ' + r.status}`);
  }
}
fs.rmSync(tmp, { recursive: true, force: true });
console.log(failed ? `\n${failed} 件失敗しました` : '\nすべて通りました');
process.exitCode = failed ? 1 : 0;
