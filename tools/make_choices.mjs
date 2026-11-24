// 選択式の選択肢を作る：入力用の答え（answer.parts[].grading）の各欄に、まちがいの選択肢を3つ足す
//   node tools/make_choices.mjs [content フォルダ]        … content/choices.json を作り直す
//   node tools/make_choices.mjs [content] --check          … 書き出さずに、今の choices.json と同じか調べる
//
// Lemma（tools/make_choices.mjs）と同じ考え方の、物理版。AI ではなく決まった手順なので、何度実行しても同じ結果。
//   候補の作り方（上から優先）：
//     0. 問題に書いてある典型的な誤答（answer.parts[].wrong。誤概念に結びついたもの）
//     1. 式の形の誤り：sin と cos の取り違え・分母と分子の逆転・2 倍／半分の付け忘れ・2乗の付け忘れ・
//        和と差の取り違え・符号（仕事・電位・変化量など、符号に意味がある量だけ）
//     2. 数の誤り：単位つきの数値は 2倍・半分・10倍（指数の1つずれ）など、有効数字の桁数をそろえて変える。
//        回・個・倍などの小さい整数は ±1 も使う。単位はそのまま残す
//   残す条件：
//     1. アプリの採点（js/answer.js の checkField）で「不正解」になること
//     2. 候補どうし・正解と、表示も値も重ならないこと
//     3. √4 や 2/4 や 1倍 のような、見ただけで不自然とわかる形でないこと
//   3つ目は、選んだ誤答をさらに1か所変えたものを優先する（正解だけが「全員の真ん中」にならないように）。
//   2つしか作れない欄は選択肢3つ、1つ以下の欄は選択肢なし（その欄は入力式）。
//
// choices.json：{ 問題ID: { spec: 入力用の答え, f: [ [正解のTeX, 誤答のTeX, …] | null, …（欄の順） ],
//                          t: [ [正解の書き方, 誤答の書き方, …] | null, … ] } }
//   spec・f は Lemma と同じ形（build_data.py は spec が今の入力用の答えと同じときだけ f を本文の ch に出す）。
//   t は検査用（tools/check_choices.mjs が、誤答がすべて不正解・正解が正解になることを採点で確かめる）。
import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { parseSpec, checkField } from '../js/answer.js';
import { freeVars } from '../js/parser.js';
import { choiceTeX, fieldInput, problemSpec, reparse, unitSuffix } from './choices_lib.mjs';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const args = process.argv.slice(2).filter(a => !a.startsWith('--'));
const content = path.resolve(args[0] || path.join(ROOT, 'content'));
const CHECK = process.argv.includes('--check');

// ---- 決まった乱数（問題ID と欄の番号から）
function seedOf(s) { let h = 2166136261; for (const c of s) { h ^= c.charCodeAt(0); h = Math.imul(h, 16777619); } return h >>> 0; }
function rng(seed) { let x = seed || 1; return () => { x ^= x << 13; x >>>= 0; x ^= x >>> 17; x ^= x << 5; x >>>= 0; return x / 4294967296; }; }

// ---- 数の書き方
const SCI = /^(-?)(\d+(?:\.\d+)?)\s*([*×])\s*10\^(?:\((-?\d+)\)|(\d+))$/;   // 6.24*10^(-11)、1.7*10^11、1.0×10^(-3)
const FRAC = /^(\()?(\d+)\/(\d+)(\))?$/;   // 9/8、(1/2)
const PLAIN_NUM = /^-?\d+(\.\d+)?$/;
function sigOf(m) {      // 有効数字の桁数（"2.0" は 2、"300" は 1、"0.25" は 2）
  const d = m.replace('-', '').replace('.', '').replace(/^0+/, '');
  return Math.max(1, m.includes('.') ? d.length : d.replace(/0+$/, '').length);
}
function fmtLike(orig, w) {
  // orig と同じ見た目（小数なら有効数字の桁数をそろえる）で w を書く。書けないときは null
  if (!(Number.isFinite(w)) || w === 0) return null;
  if (orig.includes('.')) {
    const s = w.toPrecision(sigOf(orig));
    return /e/.test(s) ? null : s;
  }
  if (Number.isInteger(w)) return String(w);
  const s = String(+w.toPrecision(Math.max(2, sigOf(orig))));
  return /e/.test(s) ? null : s;
}
function sciText(sign, mant, e, x = '*') {
  return `${sign}${mant}${x}10^${e < 0 ? `(${e})` : e}`;
}
function sciTweaks(s) {
  const m = s.match(SCI);
  if (!m) return [];
  const [, sign, mant, x, en, ep] = m;
  const e0 = Number(en ?? ep), v = Number(mant);
  const out = [];
  const add = (w, e) => {
    while (w >= 10 - 1e-9) { w /= 10; e++; }
    while (w < 1 - 1e-9) { w *= 10; e--; }
    if ((e === 0 || e === 1) && !(e0 === 0 || e0 === 1)) return;
    const t = fmtLike(mant.includes('.') || sigOf(mant) > 1 ? mant : mant, w);
    if (t) out.push(sciText(sign, t, e, x));
  };
  add(v * 2, e0); add(v / 2, e0); add(v, e0 + 1); add(v, e0 - 1); add(v * 4, e0); add(v / 4, e0); add(v * Math.SQRT2, e0);
  if (sign === '-') out.unshift(s.slice(1));
  return out;
}
const COUNT_UNIT = /^(回|個|本|枚|倍|通り|人|段|列|番目|種類)$/;
function plainNumTweaks(s, unit) {
  if (!PLAIN_NUM.test(s)) return [];
  const v = Number(s), neg = v < 0, a = Math.abs(v);
  if (a === 0) return [];
  const out = [];
  const deg = unit === '°';
  let alts;
  if (deg) alts = [a + 15, a - 15, a + 30, a - 30, a * 2, a / 2].filter(x => Number.isInteger(x) && x > 0 && x < 360);
  else if (!s.includes('.') && a <= 12 && (COUNT_UNIT.test(unit) || !unit)) alts = [a + 1, a - 1, a * 2, a % 2 === 0 ? a / 2 : a + 2, a * 3];
  else alts = [a * 2, a / 2, a * 10, a / 10, a * 4, a / 4, a * 1.5];
  for (const w of alts) {
    if (!(w > 0) || w === a) continue;
    const t = fmtLike(s.replace('-', ''), w);
    if (t && Number(t) !== a) out.push((neg ? '-' : '') + t);
  }
  if (neg) out.unshift(s.slice(1));
  return out;
}

function fracTweaks(s) {
  // 分数の答え（9/8、(1/2)倍）：分母・分子を 2倍、逆数、分子 ±1。約分できる形は looksOdd で除く
  const m = s.match(FRAC);
  if (!m || !!m[1] !== !!m[4]) return [];
  const a = +m[2], b = +m[3], par = !!m[1];
  const w = (p, q) => (q === 1 ? String(p) : par ? `(${p}/${q})` : `${p}/${q}`);
  return [w(a, b * 2), w(a * 2, b), w(b, a), w(a + 1, b), a > 1 ? w(a - 1, b) : null, w(a, b + 1), w(a * 3, b)].filter(x => x && x !== s);
}

// ---- 式の形の誤り（書式の文字列のまま変える）
function topLevel(s, chars) {
  // 括弧の外にある chars の位置（先頭の符号は除く）
  const out = []; let d = 0;
  for (let i = 0; i < s.length; i++) {
    const c = s[i];
    if (c === '(') d++;
    else if (c === ')') d--;
    else if (d === 0 && chars.includes(c) && !(c === '-' && (i === 0 || /[(^*/]/.test(s[i - 1])))) out.push(i);
  }
  return out;
}
// i の位置を囲むかっこ（無ければ全体）の中に、同じ深さの + - があるか
function inSum(s, i) {
  let a = i, d = 0;
  for (; a > 0; a--) { const c = s[a - 1]; if (c === ')') d++; else if (c === '(') { if (!d) break; d--; } }
  let b = i; d = 0;
  for (; b < s.length; b++) { const c = s[b]; if (c === '(') d++; else if (c === ')') { if (!d) break; d--; } }
  return topLevel(s.slice(a, b), '+-').length > 0;
}
const isWrapped = (x) => { if (!/^\(.*\)$/.test(x)) return false; let d = 0; for (let i = 0; i < x.length; i++) { if (x[i] === '(') d++; else if (x[i] === ')') d--; if (d === 0 && i < x.length - 1) return false; } return true; };
const unwrap = (x) => (isWrapped(x) ? x.slice(1, -1) : x);
// 分母として1かたまりのもの（(…)・文字・sqrt(…)）。μ_0 I a/(2π) log(…) のような「分数×関数」は除く
const simpleDen = (d) => isWrapped(d) || /^[A-Za-zα-ωΓΔΘΛΣΦΨΩπ](_(\{[^}]*\}|[A-Za-z0-9]+))?(\^\d+)?$/.test(d) || (/^sqrt\(/.test(d) && isWrapped(d.slice(4)));
const needsWrap = (x) => !/^[A-Za-zα-ωΓΔΘΛΣΦΨΩπ](\^\d+)?$/.test(x) && !/^[A-Za-zα-ωΔ]_(\{[^}]*\}|[A-Za-z0-9]+)$/.test(x) && !/^\d+$/.test(x) && !isWrapped(x);

function shapeTweaks(s, label) {
  const out = [];
  const sums = topLevel(s, '+-');
  // sin と cos の取り違え（1か所ずつ）
  for (const m of s.matchAll(/(?<![a-z])(sin|cos)\(/g)) {
    out.push(s.slice(0, m.index) + (m[1] === 'sin' ? 'cos' : 'sin') + s.slice(m.index + 3));
  }
  if (!sums.length && !/^-/.test(s)) {
    const divs = topLevel(s, '/');
    // 分母と分子の逆転（1つの分数のときだけ）
    if (divs.length === 1) {
      const num = s.slice(0, divs[0]), den = s.slice(divs[0] + 1);
      if (num && den && !/^\d+$/.test(num) && simpleDen(den)) {
        const n2 = unwrap(den), d2 = unwrap(num);
        out.push(`${topLevel(n2, '+-').length ? `(${n2})` : n2}/${needsWrap(d2) ? `(${d2})` : d2}`);
      }
    }
    // 2倍・半分（係数が書いてない式だけ。係数があれば数の誤りで変える）
    if (!/^\d/.test(s)) {
      out.push('2' + s);
      if (divs.length === 1) {
        const den = s.slice(divs[0] + 1), inner = unwrap(den);
        if (!/^\d/.test(inner) && simpleDen(den)) out.push(s.slice(0, divs[0] + 1) + `(2${inner})`);
      } else if (!divs.length) out.push(s + '/2');
    }
  }
  // 和・差の式は、項ごとに 2倍・半分（P_0+M g/S → P_0+2M g/S、P_0+M g/(2S)）
  if (sums.length) {
    const cuts = [0, ...sums, s.length];
    for (let k = 0; k + 1 < cuts.length; k++) {
      const a = k ? cuts[k] + 1 : 0, b = cuts[k + 1];
      const term = s.slice(a, b);
      if (!term || /^[\d-]/.test(term)) continue;
      out.push(s.slice(0, a) + '2' + term + s.slice(b));
      const dv = topLevel(term, '/');
      if (dv.length === 1) {
        const den = term.slice(dv[0] + 1);
        if (simpleDen(den) && !/^\(?\d/.test(den)) out.push(s.slice(0, a) + term.slice(0, dv[0] + 1) + `(2${unwrap(den)})` + s.slice(b));
      }
    }
  }
  // √ の中の 2 の付け忘れ（sqrt(g l) → sqrt(2g l)）
  for (const m of s.matchAll(/sqrt\((?=[A-Za-zα-ω(])/g)) out.push(s.slice(0, m.index + 5) + '2' + s.slice(m.index + 5));
  // 2乗の付け忘れ（v^2 → v）・2乗にしてしまう（1か所ずつ）
  for (const m of s.matchAll(/\^(\d+)(?![\d.])/g)) {
    const v = Number(m[1]);
    if (inSum(s, m.index)) continue;          // 和の中の項だけ次元を変えると、見ただけで誤りとわかる（a^2+y）
    if (v === 2) out.push(s.slice(0, m.index) + s.slice(m.index + m[0].length));
    else if (v >= 3) out.push(s.slice(0, m.index) + `^${v - 1}` + s.slice(m.index + m[0].length));
  }
  // 和と差の取り違え（1か所ずつ。直前が数・文字・閉じかっこのときだけ＝2項の演算）
  for (let i = 1; i < s.length; i++) {
    const c = s[i];
    if ((c === '+' || c === '-') && /[\w)πα-ωΩ]/.test(s[i - 1]) && s[i - 1] !== '^' && !/\*10\^\($/.test(s.slice(0, i))) {
      out.push(s.slice(0, i) + (c === '+' ? '-' : '+') + s.slice(i + 1));
    }
  }
  // 符号：負の答えは正にする。正の答えに - を付けるのは、符号に意味がある量（仕事・電位・エネルギーの変化など）だけ
  if (/^-/.test(s) && !topLevel(s, '+-').length) {
    out.push(s.slice(1));
    for (const x of shapeTweaks(s.slice(1), null)) if (!/^-/.test(x) && !topLevel(x, '+-').length) out.push('-' + x);
  }
  else if (!sums.length && /^(W|V|U|Δ|ΔU|ΔK|ΔE|Q|q|x|y|φ|ΔV|ΔQ|仕事|電位)/.test(label || '')) out.push('-' + s);
  return out;
}

function coefTweaks(s) {
  // 式の中の係数・分母の数（添字・指数・10^n の数は変えない）
  const out = [];
  for (const m of s.matchAll(/\d+(\.\d+)?/g)) {
    const before = s.slice(0, m.index), after = s.slice(m.index + m[0].length);
    if (/[A-Za-zα-ωΔ_{]$/.test(before) || /_\{[^}]*$/.test(before)) continue;     // 添字
    if (/\^\(?-?$/.test(before)) continue;                                          // 指数（shapeTweaks で扱う）
    if (/^\^/.test(after) && m[0] === '10') continue;                               // 10^n
    if (/[.\d]$/.test(before) || /^[.\d]/.test(after)) continue;
    const v = Number(m[0]), dec = m[0].includes('.');
    if (/\^\((\d+\/)?$/.test(before)) {                                         // 指数の分数（^(1/3)）は ±1 だけ
      for (const w of [v + 1, v - 1]) if (w > 0) out.push(before + w + after);
      continue;
    }
    const alts = dec ? [v * 2, v / 2] : [v * 2, v % 2 === 0 ? v / 2 : v + 2, ...(v <= 6 ? [v + 1, v - 1] : []), v * 3];
    for (const w of alts) {
      if (!(w > 0) || w === v) continue;
      const ws = dec ? fmtLike(m[0], w) : String(w);
      if (!ws) continue;
      let t = before + ws + after;
      // 係数の 1 は書かない（1sqrt(3) → sqrt(3)）
      if (w === 1 && /^[a-zA-Zα-ωπ(√]/.test(after) && !/[\^_/]$/.test(before)) t = before + after;
      out.push(t);
    }
  }
  if (/π/.test(s)) { out.push(s.replace(/(?<![\d.])π/, '2π')); out.push(s.replace(/(\d+)π/, (m, a) => (a === '2' ? 'π' : m))); }
  return out;
}

function rangeTweaks(s) {
  const out = [];
  for (const m of s.matchAll(/<=|>=|≦|≧|<|>/g)) {
    const swap = { '<=': '<', '>=': '>', '≦': '<', '≧': '>', '<': '<=', '>': '>=' }[m[0]];
    out.push(s.slice(0, m.index) + swap + s.slice(m.index + m[0].length));
    const rev = { '<=': '>=', '>=': '<=', '≦': '≧', '≧': '≦', '<': '>', '>': '<' }[m[0]];
    out.push(s.slice(0, m.index) + rev + s.slice(m.index + m[0].length));
  }
  return out;
}

// ---- 不自然な形を除く
function gcd(a, b) { while (b) [a, b] = [b, a % b]; return a; }
function looksOdd(s) {
  for (const m of s.matchAll(/sqrt\((\d+)\)/g)) {
    const v = +m[1]; if (v < 2) return true;
    for (let k = 2; k * k <= v; k++) if (v % (k * k) === 0) return true;   // √4, √8 のような、まだ簡単にできる形
  }
  for (const m of s.matchAll(/(?<![\w.^)])(\d+)\/(\d+)(?![\w.(])/g)) {
    const a = +m[1], b = +m[2];
    if (b === 1 || gcd(a, b) > 1) return true;
  }
  if (/(?<![\d._])0(?=[a-zA-Zα-ωπ(])|\/0(?!\.)|\^1(?!\d)|\^0(?![.\d])|\^\(1\)/.test(s)) return true;
  if (/π\/1(?!\d)|\/1(?![\d.])|log\(1\)|\(1\)\^|^1倍$|^\(1\/1\)/.test(s)) return true;
  if (/(?<![A-Za-z])\(\d+\)(?![\^\d])/.test(s.replace(/^\((\d+)\/(\d+)\)倍$/, ''))) return true;           // (2) のような余分なかっこ
  if (/[*×]10\^\(?-?[01]\)?$/.test(s)) return true;
  if (/\)\d/.test(s)) return true;
  for (const m of s.matchAll(/(?<![\w.)])(\d+)\/sqrt\((\d+)\)/g)) if (gcd(+m[1], +m[2]) > 1) return true;   // 2/sqrt(2)                                                              // sqrt(2)2π のような並び
  for (const m of s.matchAll(/(?<![\w.^_])(\d+)(?![\d.])[^/()+\-]*\/\(?(\d+)(?![\d.])/g)) if (gcd(+m[1], +m[2]) > 1) return true;   // sqrt(6G M/(2r))
  if (/\/\((π|[A-Za-zα-ω](_\w+)?)\)/.test(s)) return true;                                   // /(π) のような余分なかっこ
  const dv = topLevel(s, '/');
  if (dv.length === 1) {                                                                       // 分子と分母に同じ (…) がある
    const groups = (x) => new Set([...x.matchAll(/\(([^()]+)\)/g)].map(m => m[1]));
    const gn = groups(s.slice(0, dv[0]));
    if ([...groups(s.slice(dv[0] + 1))].some(g => gn.has(g) && /[+-]/.test(g))) return true;
  }
  // 約分できる係数（6E_1/4、2n^2h^2/(4π^2…)）
  const cf = s.match(/^-?(\d+)(?![\d.])[^/]*\/\(?(\d+)(?![\d.])/);
  if (cf && topLevel(s, '/').length === 1 && gcd(+cf[1], +cf[2]) > 1) return true;
  return false;
}

function emptyRange(g) {
  const c = g.cond;
  if (c.kind !== 'cond' || c.vars.size !== 1) return false;
  const v = [...c.vars][0];
  for (let i = -4000; i <= 4000; i++) {
    const x = i / 40;
    try { if (c.test({ [v]: x })) return false; } catch { return false; }
  }
  return true;
}
// 数値の答えどうしは 10% 以上離す
function numValue(g) {
  if (g.kind !== 'expr' || freeVars(g.value.ast).size) return null;
  const v = g.value.fn({});
  return Number.isFinite(v) ? v : null;
}
function tooClose(g, others) {
  const a = numValue(g);
  if (a === null) return false;
  return others.some(o => { const b = numValue(o.g); return b !== null && Math.abs(a - b) < 0.1 * Math.max(Math.abs(a), Math.abs(b)); });
}
function same(a, b) { const r = checkField(a, fieldInput(b)); return r.ok === true; }

// 欄 f の候補（書式の文字列。単位を除いた部分）を作る
function candidatesFor(f, wrongs) {
  const pri0 = [], pri1 = [], pri2 = [];
  for (const w of wrongs) pri0.push({ full: w });
  if (f.kind === 'expr') {
    const s = f.value.src.trim(), unit = f.value.unit || '';
    if (SCI.test(s)) pri2.push(...sciTweaks(s));
    else if (FRAC.test(s)) pri2.push(...fracTweaks(s));
    else if (PLAIN_NUM.test(s)) pri2.push(...plainNumTweaks(s, unit));
    else {
      pri1.push(...shapeTweaks(s, f.label));
      pri2.push(...coefTweaks(s));
    }
  } else if (f.kind === 'range') {
    pri1.push(...rangeTweaks(f.src)); pri2.push(...coefTweaks(f.src));
  } else if (f.kind === 'set' || f.kind === 'tuple' || f.kind === 'ratio') {
    pri2.push(...coefTweaks(f.src), ...shapeTweaks(f.src, null).filter(x => /^[{(]/.test(x) || f.kind === 'ratio'));
  }
  return { pri0, pri1, pri2 };
}

function choicesFor(f, seed, wrongs) {
  if (f.kind === 'self') return null;
  const suffix = unitSuffix(f);
  const correct = { g: f, text: f.kind === 'expr' ? f.value.src.trim() + suffix : f.src, t: choiceTeX(f) };
  const r = rng(seed);
  const { pri0, pri1, pri2 } = candidatesFor(f, wrongs);
  const mix = (xs, base) => [...new Set(xs.map(x => (typeof x === 'string' ? x : x.full)))]
    .map(x => [x, base + r()]).sort((a, b) => a[1] - b[1]).map(([x]) => x);
  // 文字列 → 欄に書く文字列（典型的な誤答は書いてあるとおり、ほかは単位を付け直す）
  const full0 = new Set(pri0.map(x => x.full));
  const cands = [...mix(pri0, 0), ...mix([...pri1, ...pri2.slice(0, 2)], 1), ...mix(pri2.slice(2), 1.5)];
  const picked = [], texs = new Set([correct.t]);
  function tryCand(c, authored) {
    // 書いてある誤答に単位が無ければ、正解と同じ単位を付ける（単位の有無で見分けられないように）
    const text = authored ? (f.kind === 'expr' && suffix && !/(\]|[^\x00-\x7F])$/.test(c) ? c + suffix : c) : (f.kind === 'expr' ? c + suffix : c);
    if (text === correct.text || (!authored && looksOdd(text))) return null;
    let g;
    try { g = reparse(f, text); } catch { return null; }
    if (!g) return null;
    const res = checkField(f, fieldInput(g));
    if (res.ok !== false || res.syntax) return null;
    if (res.error && !authored) return null;                              // 正解と区別できない・読めないものは使わない
    if (same(g, f)) return null;
    if (g.kind === 'range' && emptyRange(g)) return null;                 // 解のない範囲は使わない
    let t;
    try { t = choiceTeX(g); } catch { return null; }
    if (!t || texs.has(t)) return null;
    if (picked.some(p => same(p.g, g) || same(g, p.g))) return null;
    if (tooClose(g, [correct, ...picked])) return null;                    // 3.1 と 3.3 のような、紛らわしいだけの差は使わない
    return { g, t, text };
  }
  for (const c of cands) {
    if (picked.length >= 3) break;
    const p = tryCand(c, full0.has(c));
    if (p) { picked.push(p); texs.add(p.t); }
  }
  // 3つ目は、選んだ誤答（書いてある誤答以外）をさらに1か所変えたものを優先する
  if (picked.length >= 2 && !(picked.length === 3 && full0.has(picked[2].text))) {
    const r2 = rng(seed ^ 0x9e3779b9);
    const from = picked[Math.floor(r2() * 2)];
    if (!full0.has(from.text) && from.g.kind === f.kind) {
      const { pri1: m1, pri2: m2 } = candidatesFor(from.g, []);
      const more = [...new Set([...m1, ...m2])].map(c => [c, r2()]).sort((a, b) => a[1] - b[1]).map(x => x[0]);
      const old = picked.splice(2, 1);
      if (old.length) texs.delete(old[0].t);
      let rep = null;
      for (const c of more) { rep = tryCand(c, false); if (rep) break; }
      const back = rep || old[0];
      if (back) { picked.push(back); texs.add(back.t); }
    }
  }
  if (picked.length < 2) return null;
  return { tex: [correct.t, ...picked.map(p => p.t)], text: [correct.text, ...picked.map(p => p.text)] };
}

// ---- 実行
const pdir = path.join(content, 'problems');
const out = {};
let nProb = 0, nField = 0, nChoice = 0, nWithAll = 0, nAuto = 0;
const dirs = fs.readdirSync(pdir).sort().flatMap(u => fs.statSync(path.join(pdir, u)).isDirectory()
  ? fs.readdirSync(path.join(pdir, u)).sort().map(id => path.join(pdir, u, id, 'problem.json')).filter(p => fs.existsSync(p)) : []);
for (const file of dirs) {
  const p = JSON.parse(fs.readFileSync(file, 'utf8'));
  if (p.quality?.status !== 'published') continue;
  const spec = problemSpec(p);
  if (!spec) continue;
  let sp;
  try { sp = parseSpec(spec); } catch { continue; }
  nAuto++;
  const fields = [], wrongs = [];
  // 典型的な誤答（answer.parts[].wrong）は、欄が1つの小問だけで使う
  const partWrong = Object.fromEntries(p.answer.parts.filter(x => x.grading).map(x => [x.id === 'main' ? '' : x.id, (x.wrong || []).map(w => w.value)]));
  for (const part of sp.parts) for (const fl of part.fields) {
    fields.push(fl);
    wrongs.push(part.fields.length === 1 ? partWrong[part.n] || [] : []);
  }
  const res = fields.map((fl, k) => { if (fl.kind !== 'self') nField++; return choicesFor(fl, seedOf(p.id + '#' + k), wrongs[k]); });
  const got = res.filter(Boolean).length;
  if (!got) continue;
  nProb++; nChoice += got;
  if (got === fields.filter(x => x.kind !== 'self').length) nWithAll++;
  out[p.id] = { spec, f: res.map(x => (x ? x.tex : null)), t: res.map(x => (x ? x.text : null)) };
}
const ids = Object.keys(out).sort();
const text = '{\n' + ids.map(id => JSON.stringify(id) + ':' + JSON.stringify(out[id])).join(',\n') + '\n}\n';
const dest = path.join(content, 'choices.json');
if (CHECK) {
  const ok = fs.existsSync(dest) && fs.readFileSync(dest, 'utf8') === text;
  console.log('content/choices.json：', ok ? '同じ' : '違いがあります → node tools/make_choices.mjs');
  process.exit(ok ? 0 : 1);
}
fs.writeFileSync(dest, text);
console.log(`選択肢：${nProb}問（入力して採点できる ${nAuto}問・自動採点の欄 ${nField} のうち ${nChoice} 欄、すべての欄に選択肢がある問題 ${nWithAll}問） → ${path.relative(process.cwd(), dest)}`);
