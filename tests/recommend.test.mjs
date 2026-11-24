// 「今日の演習」（catalog.js の pickToday）の推薦の性質を確かめる（v0.6。Lemma v0.21 の tests/recommend.test.mjs から解法カードの部分を除いたもの）
//   node tests/recommend.test.mjs
// 前半：架空の学習者で40日間「今日の演習」を回し、推薦の規則が守られているか（学習成果ではなく規則だけを見る。正答率は仮定）
// 後半：境界条件（問題数・難易度の選択・復習が大量・選べる問題が少ない）で、制限を越えず、足りないときは理由を返すか
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
globalThis.fetch = async (url) => {
  const p = path.join(root, String(url));
  if (!fs.existsSync(p)) return { ok: false, status: 404, json: async () => { throw new Error('404'); } };
  const text = fs.readFileSync(p, 'utf8');
  return { ok: true, status: 200, json: async () => JSON.parse(text), text: async () => text };
};
const CAT = await import('../js/catalog.js');
const { review } = await import('../js/scheduler.js');
await CAT.loadCatalog();

let pass = 0, fail = 0;
const ok = (c, m) => { if (c) pass++; else { fail++; console.log('  FAIL ' + m); } };
const DAY = 86400000;
const RULES = { due: /^(前回完答できなかった問題の復習|復習の時期)$/, weak: /^苦手な単元「.+」を補う$/, fresh: /^(新しい問題|解いてから日がたった問題（選んだ難易度の未着手が残っていない）)$/ };

// ---- 架空の学習者（正答率は仮定）
const PERSONAS = {
  steady: { name: '毎日8問・標準まで', levels: [0, 1, 2], size: 8, p: (e) => 0.75 - 0.1 * e.lv },
  strong: { name: '入試レベル・正答率高め', levels: [2, 3, 4], size: 6, p: () => 0.85 },
  weakMec: { name: '力学が苦手', levels: [0, 1, 2, 3], size: 8, p: (e) => (e.f.startsWith('MEC') ? 0.2 : 0.8) },
  small: { name: '1日4問・基礎だけ', levels: [0], size: 4, p: () => 0.5 },
};
function seeded(seed) { let x = seed >>> 0; return () => { x = (x * 1664525 + 1013904223) >>> 0; return x / 4294967296; }; }
function weakest(st) {
  let w = null;
  for (const f of CAT.FIELDS) {
    const last = new Map(); for (const l of st.log) { const e = CAT.entryById.get(l.p); if (e && e.f === f.key) last.set(l.p, { l, e }); }
    if (last.size < 5) continue;
    let sw = 0, ss = 0; for (const { l, e } of last.values()) { const wt = 1 + 0.25 * e.lv; sw += wt; ss += wt * (l.g >= 3 ? 1 : l.g === 2 ? 0.5 : 0); }
    const sc = (ss + 1) / (sw + 2);
    if (sc < 0.5 && (!w || sc < w.score)) w = { f, score: sc };
  }
  return w;
}
const realNow = Date.now;
for (const [k, P] of Object.entries(PERSONAS)) {
  const rnd = seeded(k.length * 31 + 2);
  const st = { problems: {}, pstat: {}, log: [] };
  const t0 = Date.UTC(2026, 0, 5, 9);
  let minFresh = Infinity, levelViolations = 0, repeats3d = 0, fieldMaxShare = 0, mismatches = [], days = 0;
  for (let d = 0; d < 40; d++) {
    const now = t0 + d * DAY;
    Date.now = () => now;
    const w = weakest(st);
    const plan = CAT.pickToday(st, { levels: P.levels, size: P.size, weak: w ? w.f.key : null, weakName: w ? w.f.name : '' });
    if (!plan.steps.length) continue;
    days++;
    const nonDue = plan.items.filter(x => x.tag !== 'due');
    minFresh = Math.min(minFresh, nonDue.length);
    for (const x of plan.items) {
      const e = CAT.entryById.get(x.pid);
      if (x.tag !== 'due' && !P.levels.includes(e.lv)) levelViolations++;
      if (x.tag !== 'due' && st.problems[x.pid] && now - st.problems[x.pid].last < 3 * DAY) repeats3d++;
      if (!RULES[x.tag] || !RULES[x.tag].test(x.why)) mismatches.push(`${x.tag}:${x.why}`);
      if (x.tag === 'due' && !(CAT.dueEntries(st, {}, now).some(y => y.id === x.pid))) mismatches.push(`due だが復習の時期でない ${x.pid}`);
    }
    const byF = {}; for (const x of plan.items) { const f = CAT.entryById.get(x.pid).f; byF[f] = (byF[f] || 0) + 1; }
    if (plan.items.length >= 4) fieldMaxShare = Math.max(fieldMaxShare, Math.max(...Object.values(byF)) / plan.items.length);
    // 解く
    plan.steps.forEach((s, i) => {
      const e = CAT.entryById.get(s.pid), t = now + i * 600000;
      const g = rnd() < P.p(e) ? 3 : rnd() < 0.5 ? 2 : 1;
      st.pstat[s.pid] = review(st.pstat[s.pid], g, t);
      st.problems[s.pid] = { last: t, n: ((st.problems[s.pid] || {}).n || 0) + 1 };
      st.log.push({ p: s.pid, t, g });
    });
  }
  Date.now = realNow;
  const n = P.name;
  ok(days >= 39, `${n}：40日とも出題できる（${days}日）`);
  ok(minFresh >= Math.min(2, P.size), `${n}：毎回、復習以外の問題が2問以上（最少 ${minFresh}）`);
  ok(levelViolations === 0, `${n}：選んだ難易度の外から新しい問題を出さない（${levelViolations}）`);
  ok(repeats3d === 0, `${n}：復習以外で3日以内に同じ問題を出さない（${repeats3d}）`);
  ok(fieldMaxShare <= 0.5, `${n}：1回で同じ単元が半分を超えない（${Math.round(fieldMaxShare * 100)}%）`);
  ok(mismatches.length === 0, `${n}：推薦の理由と選んだ規則が一致（${mismatches.slice(0, 2)}）`);
}

// ---- 境界条件（Lemma v0.21）：制限を越えず、足りないときは理由を返す
{
  const now = Date.now(), old = now - 90 * DAY;
  const fresh = () => ({ problems: {}, pstat: {}, log: [] });
  const withDue = (es) => { const s = fresh(); for (const e of es) { s.problems[e.id] = { last: old, n: 1 }; s.pstat[e.id] = { S: 1, D: 5, last: old, reps: 1, lapses: 0 }; s.log.push({ p: e.id, t: old, g: 3 }); } return s; };
  const check = (plan, size, levels, label) => {
    const ids = plan.items.map(x => x.pid);
    ok(new Set(ids).size === ids.length, `${label}：同じ問題を2回入れない`);
    ok(plan.items.every(x => x.tag === 'due' || !levels || levels.includes(CAT.entryById.get(x.pid).lv)), `${label}：選んだ難易度を越えない（復習を除く）`);
    const nd = {}, all = {};
    for (const x of plan.items) { const f = CAT.entryById.get(x.pid).f; all[f] = (all[f] || 0) + 1; if (x.tag !== 'due') nd[f] = (nd[f] || 0) + 1; }
    ok(Math.max(0, ...Object.values(nd)) <= Math.max(2, Math.ceil(size / 3)) && Math.max(0, ...Object.values(all)) <= Math.max(2, Math.ceil(size / 2)), `${label}：単元の上限を越えない`);
    ok(plan.steps.length === size || (plan.short && plan.short.got === plan.steps.length && plan.short.text), `${label}：${size}問に足りなければ理由を返す（${plan.steps.length}問${plan.short ? '・' + plan.short.cause : ''}）`);
  };
  for (const size of [1, 3, 5, 12, 20]) check(CAT.pickToday(fresh(), { size, levels: [0, 1, 2, 3, 4] }), size, [0, 1, 2, 3, 4], `size ${size}`);
  const p0 = CAT.pickToday(fresh(), { size: 8 });
  ok(p0.count.due === undefined && p0.steps.length === 8, '復習待ちが無いときは新しい問題だけで8問');
  const many = CAT.ENTRIES.filter((_, i) => i % 3 === 0).slice(0, 200);
  const pm = CAT.pickToday(withDue(many), { size: 8 });
  check(pm, 8, null, '復習待ち200問');
  ok(pm.count.due === 6, `復習待ちが大量でも新しい学習の枠を2問残す（復習 ${pm.count.due}）`);
  const few = CAT.pickToday(withDue(many.slice(0, 5)), { size: 8 });
  ok(few.count.due === 5, '復習が少ないときは新しい問題の枠を3問残す（復習5問はすべて出す）');
  check(CAT.pickToday(fresh(), { size: 8, levels: [5] }), 8, [5], '物理探究だけ');
  const lv0 = CAT.ENTRIES.filter(e => e.lv === 0);
  const almost = fresh(); for (const e of lv0) { almost.problems[e.id] = { last: now - 3600000, n: 1 }; almost.pstat[e.id] = { S: 30, D: 5, last: now - 3600000, reps: 2, lapses: 0 }; }
  const pa = CAT.pickToday(almost, { size: 8, levels: [0] });
  check(pa, 8, [0], '基礎をすべて解き終えた直後');
  ok(pa.steps.length === 0 && pa.short && pa.short.cause === 'pool', '基礎をすべて今日解いたなら、出せる問題が無い理由を返す');
  const f0 = CAT.ENTRIES[0].f, one = fresh();
  for (const e of CAT.ENTRIES.filter(e => e.f !== f0 && e.lv === 0)) { one.problems[e.id] = { last: now - 60000, n: 1 }; one.pstat[e.id] = { S: 30, D: 5, last: now - 60000, reps: 2, lapses: 0 }; }
  const p1 = CAT.pickToday(one, { size: 8, levels: [0] });
  check(p1, 8, [0], '未着手が1単元だけ');
  ok(p1.short && p1.short.cause === 'field', '未着手が1単元だけなら、偏らせないために少なくした理由を返す');
}

console.log(`recommend.test: ${pass} 件 OK${fail ? `、${fail} 件失敗` : ''}`);
process.exit(fail ? 1 : 0);
