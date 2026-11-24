// 選択式の選択肢（content/choices.json）を、アプリの採点（js/answer.js）で確かめる
//   node tools/check_choices.mjs [content フォルダ]
// 確かめること（問題ごと・欄ごと）：
//   - spec が今の入力用の答えと同じ（違えば古い選択肢。node tools/make_choices.mjs で作り直す）
//   - f（表示する TeX）と t（書き方）の長さが欄の数と同じで、null の位置もそろっている
//   - 正解（t の先頭）を入れると正解になる。表示（f の先頭）が欄の正解の TeX と同じ
//   - 誤答（t の2つ目以降）はすべて不正解になり、表示が f と一致する
//   - 選択肢どうしの表示が重ならず、値も互いに等しくない
//   - 選択肢は 3つか4つ（誤答 2〜3）
//   - KaTeX で描画できる
import fs from 'fs';
import path from 'path';
import { createRequire } from 'module';
import { fileURLToPath } from 'url';
import { parseSpec, checkField } from '../js/answer.js';
import { choiceTeX, fieldInput, problemSpec, reparse } from './choices_lib.mjs';

const require = createRequire(import.meta.url);
const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const content = path.resolve(process.argv[2] || path.join(ROOT, 'content'));
const katex = require(path.join(ROOT, 'vendor', 'katex.min.js'));

const choices = JSON.parse(fs.readFileSync(path.join(content, 'choices.json'), 'utf8'));
const specs = {};
const pdir = path.join(content, 'problems');
for (const u of fs.readdirSync(pdir)) {
  if (!fs.statSync(path.join(pdir, u)).isDirectory()) continue;
  for (const id of fs.readdirSync(path.join(pdir, u))) {
    const f = path.join(pdir, u, id, 'problem.json');
    if (!fs.existsSync(f)) continue;
    const p = JSON.parse(fs.readFileSync(f, 'utf8'));
    if (p.quality?.status === 'published') specs[p.id] = problemSpec(p);
  }
}

const bad = [];
let nProb = 0, nField = 0, nOpt = 0, stale = 0;
for (const [id, c] of Object.entries(choices)) {
  if (!(id in specs)) { bad.push(`${id}: 問題がない（公開されていない）`); continue; }
  if (specs[id] !== c.spec) { stale++; bad.push(`${id}: 入力用の答えが変わっている（古い選択肢）`); continue; }
  let sp;
  try { sp = parseSpec(c.spec); } catch (e) { bad.push(`${id}: 入力用の答えが読めない ${e.message}`); continue; }
  const fields = sp.parts.flatMap(p => p.fields);
  if (!Array.isArray(c.f) || !Array.isArray(c.t) || c.f.length !== fields.length || c.t.length !== fields.length) {
    bad.push(`${id}: f・t の長さが欄の数（${fields.length}）と違う`); continue;
  }
  nProb++;
  fields.forEach((f, k) => {
    const F = c.f[k], T = c.t[k], at = `${id}#${k + 1}`;
    if (F === null || T === null) { if (F !== T) bad.push(`${at}: f と t の null がそろっていない`); return; }
    if (f.kind === 'self') { bad.push(`${at}: 自己採点の欄に選択肢がある`); return; }
    nField++;
    if (F.length !== T.length || F.length < 3 || F.length > 4) { bad.push(`${at}: 選択肢の数が正しくない（${F.length}）`); return; }
    nOpt += F.length;
    if (F[0] !== choiceTeX(f)) bad.push(`${at}: 正解の表示が欄の正解と違う ${F[0]}`);
    const gs = [];
    T.forEach((text, i) => {
      let g;
      try { g = reparse(f, text); } catch (e) { g = null; }
      if (!g) { bad.push(`${at}: 「${text}」が読めない`); return; }
      const r = checkField(f, fieldInput(g));
      if (i === 0 && r.ok !== true) bad.push(`${at}: 正解「${text}」が正解にならない`);
      if (i > 0 && r.ok !== false) bad.push(`${at}: 誤答「${text}」が不正解にならない（${JSON.stringify(r)}）`);
      if (i > 0 && r.syntax) bad.push(`${at}: 誤答「${text}」が書式の誤りになる`);
      if (choiceTeX(g) !== F[i]) bad.push(`${at}: 「${text}」の表示が f と違う ${F[i]}`);
      try { katex.renderToString(F[i], { throwOnError: true, displayMode: true, strict: false }); }
      catch (e) { bad.push(`${at}: KaTeX で描けない ${F[i]} ${e.message.slice(0, 80)}`); }
      gs.push(g);
    });
    if (new Set(F).size !== F.length) bad.push(`${at}: 表示が重なる選択肢がある`);
    for (let i = 0; i < gs.length; i++) for (let j = i + 1; j < gs.length; j++) {
      if (checkField(gs[i], fieldInput(gs[j])).ok === true) bad.push(`${at}: 「${T[i]}」と「${T[j]}」が同じ値`);
    }
  });
}
// 選択肢があるはずなのに choices.json に無い問題（作り直し忘れ）
const missing = Object.entries(specs).filter(([id, s]) => s && !(id in choices)).map(([id]) => id);
console.log(`選択肢：${nProb}問・${nField}欄・${nOpt}個を採点で確かめました（入力して採点できる問題 ${Object.values(specs).filter(Boolean).length}問、選択肢のない問題 ${missing.length}問）`);
for (const b of bad.slice(0, 40)) console.log('  ✗', b);
if (bad.length > 40) console.log(`  …ほか ${bad.length - 40}件`);
if (stale) console.log(`  → 入力用の答えが変わった問題が ${stale}問あります。node tools/make_choices.mjs で作り直してください`);
process.exitCode = bad.length ? 1 : 0;
