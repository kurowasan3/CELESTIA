// 選択式の選択肢の共通処理（tools/make_choices.mjs と tools/check_choices.mjs が使う）
import { parseSpec, fieldTeX, exprTeX } from '../js/answer.js';

// 問題（content/problems/*/*/problem.json）の入力用の答え。build_data.py の grading_spec() と同じ
export function problemSpec(p) {
  const lines = [];
  for (const x of p.answer?.parts || []) {
    if (!('grading' in x)) continue;
    lines.push(x.id === 'main' ? x.grading : `(${x.id}) ${x.grading}`);
  }
  return lines.length ? lines.join('\n') : null;
}

export function lineFor(f, text) { return (f.label ? `${f.label}：` : '') + text; }

// 1つの欄として読み直す（種類が変わったら null）
export function reparse(f, text) {
  const sp = parseSpec(lineFor(f, text));
  if (!sp || sp.parts.length !== 1 || sp.parts[0].fields.length !== 1) return null;
  const g = sp.parts[0].fields[0];
  return g.kind === f.kind ? g : null;
}

// アプリの入力欄に打ち込む文字列（単位は欄の横に出るので打たない）
export function fieldInput(g) {
  if (g.kind === 'cases') return { mode: 'cases', rows: g.cases.map(c => ({ cond: c.condSrc, val: c.valSrc })) };
  if (g.kind === 'set') return { mode: 'plain', text: g.items.map(i => i.src).join(', ') };
  if (g.kind === 'expr') return { mode: 'plain', text: g.value.src };
  return { mode: 'plain', text: g.src };
}

// 欄の書き方のうち、値のあとの単位の部分（" [m/s]"・"倍" など）
export function unitSuffix(f) {
  if (f.kind !== 'expr') return '';
  const src = f.src.trim(), v = f.value.src.trim();
  return src.startsWith(v) ? src.slice(v.length) : (f.value.unit ? ` [${f.value.unit}]` : '');
}

// ---- 選択肢の TeX
// 数値は書いたとおりの桁数で（2.0 を 2 にしない）、10 の累乗は ×10^n で、単位は立体で
const SCI = /^(-?)(\d+(?:\.\d+)?)\s*[*×]\s*10\^(?:\((-?\d+)\)|(\d+))$/;
export function numTeX(s) {
  s = s.trim();
  const m = s.match(SCI);
  if (m) return `${m[1]}${m[2]}\\times 10^{${m[3] ?? m[4]}}`;
  if (/^-?\d+(\.\d+)?$/.test(s)) return s;
  return null;
}
const JP = /[^\x00-\x7FΩμ·°℃α-ωΑ-Ω]/;
export function unitTeX(u) {
  if (!u) return '';
  if (u === '°') return '^\\circ';
  if (u === '℃') return '\\,{}^\\circ\\mathrm{C}';
  if (JP.test(u)) return `\\text{${u}}`;
  const t = u.replace(/\^\(([^)]*)\)/g, '^{$1}').replace(/\^(-?\d+)/g, '^{$1}')
    .replace(/Ω/g, '\\Omega ').replace(/μ/g, '\\mu ').replace(/·/g, '\\cdot ').replace(/℃/g, '{}^\\circ C').replace(/ +$/, '');
  return `\\,\\mathrm{${t}}`;
}
export function choiceTeX(g) {
  if (g.kind === 'expr') {
    const s = g.value.src.trim();
    return spaceCommands((numTeX(s) ?? exprTeX(s) ?? s) + unitTeX(g.value.unit));
  }
  return spaceCommands(fieldTeX(g));
}
// toTeX は \pi のすぐ後ろに文字を並べることがある（4\pik_{0}）。命令名と文字の間に空白を入れる
const CMDS = 'alpha|beta|gamma|delta|epsilon|varepsilon|zeta|eta|theta|iota|kappa|lambda|mu|nu|xi|pi|rho|sigma|tau|upsilon|phi|varphi|chi|psi|omega|Gamma|Delta|Theta|Lambda|Xi|Pi|Sigma|Phi|Psi|Omega|cdot|times';
export function spaceCommands(t) {
  return t.replace(new RegExp(`\\\\(${CMDS})(?=[A-Za-z])`, 'g'), '\\$1 ');
}
