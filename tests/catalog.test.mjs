// catalog.js の単体テスト（ブラウザなしで data/ を読む）
//   node tests/catalog.test.mjs
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
// catalog.js は相対 URL で fetch するので、リポジトリのファイルを返す fetch に差し替える
globalThis.fetch = async (url) => {
  const p = path.join(root, String(url));
  if (!fs.existsSync(p)) return { ok: false, status: 404, json: async () => { throw new Error('404'); } };
  const text = fs.readFileSync(p, 'utf8');
  return { ok: true, status: 200, json: async () => JSON.parse(text), text: async () => text };
};
const CAT = await import('../js/catalog.js');

let fail = 0, pass = 0;
const ok = (cond, msg) => { if (cond) pass++; else { fail++; console.log('FAIL', msg); } };

const meta = await CAT.loadCatalog();
const idx = JSON.parse(fs.readFileSync(path.join(root, 'data/index.json'), 'utf8'));
ok(meta.ok && meta.format === 'celestia.data/1', '問題データを読めること');
ok(CAT.ENTRIES.length === idx.count, `ENTRIES の数 ${CAT.ENTRIES.length} = ${idx.count}`);
ok(meta.warnings.length === 0, `形式の警告がないこと ${meta.warnings}`);
ok(CAT.FIELDS.every(f => f.subject && f.fieldKey && /^[A-Z]{3}-\d{2}$/.test(f.key)), '単元は「分野 → 単元」の形');
ok(CAT.LEVELS[2] === '入試基礎' && CAT.SUBJECT_ORDER[0] === '力学', '難易度と分野の並びをデータから読む');
const tax = JSON.parse(fs.readFileSync(path.join(root, 'content/taxonomy.json'), 'utf8'));
const emptyUnits = tax.units.filter(u => !CAT.ENTRIES.some(e => e.f === u.key)).map(u => u.key);
ok(emptyUnits.length === 0, `すべての単元に問題がある（ない単元：${emptyUnits}）`);
const fieldsWithSim = new Set(JSON.parse(fs.readFileSync(path.join(root, 'content/simulations.json'), 'utf8')).simulations.map(s => s.field));
ok(tax.fields.filter(f => !['MSY', 'ESY'].includes(f.key)).every(f => fieldsWithSim.has(f.key)), 'すべての分野（総合を除く）にシミュレーションがある');

ok(CAT.pool({ levels: [0] }).length > 0, 'pool() で基礎の問題がとれる');
const st = JSON.parse(fs.readFileSync(path.join(root, 'tests/fixtures/state-v1.json'), 'utf8'));
const st05 = JSON.parse(fs.readFileSync(path.join(root, 'tests/fixtures/state-v0.5.json'), 'utf8'));
ok([...st05.log.map(l => l.p), ...Object.keys(st05.pstat)].every(id => CAT.entryById.has(id)), 'v0.5 の記録の fixture の問題ID がすべて今のデータにある');
ok(CAT.pickPractice(st, { size: 3 }).length === 3, 'pickPractice()');
ok(Array.isArray(CAT.pickReview(st, { size: 5 })), 'pickReview()');

// 今日の演習（v0.6。Lemma v0.16〜v0.21 の推薦から解法カードを除いたもの）：復習待ちを先に、残りを新しい問題で。苦手な単元から最大2問
{
  const now = Date.now(), old = now - 60 * 86400000;
  const ids = CAT.ENTRIES.slice(0, 3).map(e => e.id);
  const s1 = { problems: {}, pstat: {}, log: [] };
  for (const id of ids) { s1.problems[id] = { last: old, n: 1 }; s1.pstat[id] = { S: 1, D: 5, last: old, reps: 1, lapses: 0 }; }
  const weak = CAT.FIELDS.find(f => f.key !== CAT.entryById.get(ids[0]).f).key;
  const plan = CAT.pickToday(s1, { size: 8, weak, weakName: 'x' });
  const got = plan.steps.map(x => x.pid);
  ok(plan.due === 3 && plan.fresh === 5 && got.length === 8 && new Set(got).size === 8, `今日の演習：復習3問＋新しい問題5問（${plan.due}+${plan.fresh}）`);
  ok(ids.every(id => got.includes(id)), '今日の演習に復習待ちがすべて入る');
  ok(got.filter(id => CAT.entryById.get(id).f === weak && !ids.includes(id)).length >= 2, '今日の演習に苦手な単元の新しい問題が2問入る');
  const RULES = [/^前回完答できなかった問題の復習$/, /^復習の時期$/, /^苦手な単元「.+」を補う$/, /^新しい問題$/, /^解いてから日がたった問題（選んだ難易度の未着手が残っていない）$/];
  ok(plan.steps.every(x => RULES.some(r => r.test(x.why))), '今日の演習：すべての問題に規則どおりの理由');
  const only = CAT.pickToday({ problems: {}, pstat: {}, log: [] }, { size: 6, levels: [0] });
  ok(only.due === 0 && only.steps.length === 6 && only.steps.every(x => CAT.entryById.get(x.pid).lv === 0), '今日の演習：難易度の選択を守る');
  ok(CAT.pickToday(st, { size: 8 }).steps.length === 8, '今日の演習：記録のある端末（fixture）でも出題できる');
}
// 類題：同じ単元・難易度の差1以内・共通のタグ・概念・法則あり
{
  const e0 = CAT.ENTRIES.find(e => (e.tg || []).length >= 2 && (e.c || []).length);
  const sim = CAT.similarTo({ problems: {} }, e0.id);
  const share = (a, b) => (a || []).some(t => (b || []).includes(t));
  ok(sim && sim.id !== e0.id && sim.f === e0.f && Math.abs(sim.lv - e0.lv) <= 1 && (share(sim.tg, e0.tg) || share(sim.c, e0.c) || share(sim.l, e0.l) || share(sim.m, e0.m)), `類題：${e0.id} → ${sim && sim.id}`);
  const sim2 = CAT.similarTo({ problems: {} }, e0.id, [sim.id]);
  ok(!sim2 || sim2.id !== sim.id, '類題：除外した問題は選ばない');
  ok(CAT.similarTo({ problems: {} }, 'no-such-id') === null, '類題：無い問題は null');
}
// 頻出度（v0.6）：index.json の v は 1・3 だけ（2 は既定で書かない）。無い問題は中（2）。データに頻出度があれば、高い問題が多めに出る
ok(CAT.ENTRIES.every(e => e.v === undefined || e.v === 1 || e.v === 3), '頻出度は 1 か 3 だけが書かれる（2 は既定）');
ok(CAT.valueOf({}) === 2 && CAT.valueOf({ v: 3 }) === 3 && CAT.valueOf({ v: 9 }) === 2 && CAT.VALUES[3] === '高', 'valueOf：無い・範囲外は 2（中）');
if (CAT.ENTRIES.some(e => e.v === 3) && CAT.ENTRIES.some(e => e.v === 1)) {
  const lv = [...new Set(CAT.ENTRIES.map(e => e.lv))].sort((a, b) => CAT.pool({ levels: [b] }).length - CAT.pool({ levels: [a] }).length)[0];
  const hi = CAT.pool({ levels: [lv] }).filter(e => CAT.valueOf(e) === 3).length / CAT.pool({ levels: [lv] }).length;
  let n = 0, h = 0;
  for (let i = 0; i < 40; i++) for (const it of CAT.pickPractice({ problems: {} }, { size: 10, levels: [lv] })) { n++; if (CAT.valueOf(CAT.entryById.get(it.pid)) === 3) h++; }
  ok(h / n > hi, `頻出度の高い問題が多めに出る（出題 ${(h / n).toFixed(2)} > 全体 ${hi.toFixed(2)}）`);
} else {
  // data に頻出度がまだ無いとき：一時的に3分の1の問題を高（3）にして、重みつきの選び方だけを確かめる（あとで元に戻す）
  const pool0 = CAT.pool({ levels: [0] });
  pool0.forEach((e, i) => { if (i % 3 === 0) e.v = 3; });
  let n = 0, h = 0;
  for (let i = 0; i < 200; i++) for (const it of CAT.pickPractice({ problems: {} }, { size: 10, levels: [0] })) { n++; if (CAT.valueOf(CAT.entryById.get(it.pid)) === 3) h++; }
  pool0.forEach(e => { delete e.v; });
  ok(h / n > 0.4, `頻出度の高い問題が多めに出る（全体の3分の1 → 出題 ${(h / n).toFixed(2)}。data に v が無いので一時的に付けて確認）`);
}

// 全問の本文を読み、共通ビューが作れること
await CAT.ensureBodies(CAT.ENTRIES.map(e => e.id));
ok(CAT.staleDataFiles().length === 0, `本文ファイルの版が揃っていること ${CAT.staleDataFiles()}`);
let bad = 0;
for (const e of CAT.ENTRIES) {
  const p = CAT.getProblem(e.id);
  const good = p && p.id === e.id && p.field === e.f && p.level === e.lv && p.subject && p.fieldName && p.html && p.html.q
    && (!!p.grading === !!e.in) && (!p.grading || typeof p.grading.spec === 'string') && p.figs && typeof p.figs === 'object';
  if (!good) { bad++; if (bad < 5) console.log('  bad', e.id, p && Object.keys(p)); }
}
ok(bad === 0, `getProblem() が全 ${CAT.ENTRIES.length}問で正しい形（だめ ${bad}問）`);
ok(CAT.getProblem('no-such-id') === undefined, '無い ID は undefined');
// 選択式の選択肢（v0.6）：index.json の ch: 1 と本文の ch（入力欄と同じ並び・各欄は [正解, 誤答…] か null）が揃っている
{
  const { parseSpec } = await import('../js/answer.js');
  let chBad = 0, chN = 0;
  for (const e of CAT.ENTRIES) {
    const b = CAT.getBody(e.id);
    if (!!e.ch !== Array.isArray(b.ch)) { chBad++; if (chBad < 5) console.log('  ch の印が違う', e.id); continue; }
    if (!b.ch) continue;
    chN++;
    const nf = parseSpec(b.in).parts.reduce((a, p) => a + p.fields.length, 0);
    if (b.ch.length !== nf || !b.ch.some(Boolean) || !b.ch.every(x => x === null || (Array.isArray(x) && x.length >= 2 && x.every(t => typeof t === 'string')))) { chBad++; if (chBad < 5) console.log('  ch の形', e.id); }
    if (CAT.getProblem(e.id).ch !== b.ch) chBad++;
  }
  ok(chBad === 0, `選択肢のある問題 ${chN}問：一覧の印と本文の形が揃う（だめ ${chBad}）`);
}

ok(CAT.getBody('_meta') === undefined, '_meta を本文として登録しないこと');

// 図：本文の枠と figs の対応。問題文に解答用の図の枠が無いこと
let figBad = 0;
for (const e of CAT.ENTRIES) {
  const b = CAT.getBody(e.id), figs = b.figs || {};
  const refs = (html) => [...html.matchAll(/data-fig="([^"]+)"/g)].map(m => m[1]);
  for (const r of refs(b.q)) if (!figs[r] || figs[r].role === 'answer') { figBad++; console.log('  問題文の図', e.id, r); }
  for (const r of refs(b.s)) if (!figs[r]) { figBad++; console.log('  解説の図', e.id, r); }
  for (const f of Object.values(figs)) if (f.src && !fs.existsSync(path.join(root, f.src))) { figBad++; console.log('  ファイル', e.id, f.src); }
}
ok(figBad === 0, '図の枠・情報・ファイルがそろい、問題文に解答用の図が無い');

// 概念・法則・誤概念
const K = await CAT.loadKnowledge();
ok(K.concepts.length && K.concept['normal-force'] && K.law['newton-2nd'] && K.misconception['MC-NORMAL-EQ-MG'], 'knowledge.json を読める');
ok(CAT.ENTRIES.every(e => (e.c || []).every(c => K.concept[c]) && (e.m || []).every(m => K.misconception[m])), '問題の概念・誤概念がすべて knowledge.json にある');

// 学習記録のサンプルが今のデータで使えること
const ids = new Set([...Object.keys(st.pstat), ...Object.keys(st.problems), ...st.log.map(l => l.p)]);
const missing = [...ids].filter(i => !CAT.entryById.has(i));
ok(missing.length === 0, `記録にある問題ID がすべてある（無い ${missing.slice(0, 5)}）`);
ok(st.practice.fields.every(k => CAT.fieldByKey[k]), '記録にある単元ID がすべてある');

// 検索・問題セット
ok(CAT.normText('最大・最小') === CAT.normText('最大最小') && CAT.normText('ＡＢ　c') === 'abc', 'normText が表記ゆれを吸収する');
const r1 = await CAT.searchProblems('斜面');
ok(r1.length >= 2 && r1.some(e => e.id === 'MEC-INC-001'), `「斜面」で見つかる（${r1.length}問）`);
const r2 = await CAT.searchProblems('垂直抗力');
ok(r2.some(e => e.id === 'MEC-CON-001'), '概念名「垂直抗力」でも見つかる');
const r3 = await CAT.searchProblems('オームの法則', { levels: [0] });
ok(r3.length > 0 && r3.every(e => e.lv === 0), `法則名・難易度の絞り込み（${r3.length}問）`);
ok((await CAT.searchProblems('')).length === 0, '空の検索は0件');
const cols = await CAT.loadCollections();
ok(cols.length >= 1 && cols.every(c => c.items.every(it => CAT.entryById.has(it.id))), `問題セット ${cols.length}件、すべての問題がある`);
const pc = CAT.pickCollection({ problems: { [cols[0].items[0].id]: { last: 1, n: 1 } } }, cols[0], 3);
ok(pc.length === 3 && pc[0].pid === cols[0].items[1].id, 'pickCollection は未着手の問題から順に');

console.log(`catalog.test: ${pass} 件 OK${fail ? `、${fail} 件失敗` : ''}`);
process.exitCode = fail ? 1 : 0;
