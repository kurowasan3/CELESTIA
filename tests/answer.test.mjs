import { parseSpec, checkField, fieldTeX } from '../js/answer.js';
// 自動採点（answer.js）の回帰テスト
let passed = 0;
const spec = parseSpec(`(1) 14/3
(2) 最大値：7/4 ; 最小値：-3-2sqrt(3)
(3) {1, 35, 251, 505}
(4) (30/11, 48/11)
(5) 9:4
(6) a<-3 または 3<a
(7) 場合分け{ m<=1/2 → 6-2m | m>=1/2 → (14m+8)/3 }
(8) 証明
(9) 18個
(10) -16<k<9
(11) {4+sqrt(5), 4-sqrt(5)}
(12) 120°
(13) $y$：9x-16
(14) 3*2^(n-1)-n-1
(15) なし
(16) 場合分け{ a<0 → a | 0<=a<=2 → -a^2+a | a>2 → 4-3a }`);
const P = (n) => spec.parts.find(p => p.n == n).fields;
const T = (n, i, inp, want) => { const r = checkField(P(n)[i], typeof inp === 'string' ? { mode: 'plain', text: inp } : inp); const pass = r.ok === want; if (!pass) console.log('FAIL', n, JSON.stringify(inp), r); else passed++; if(!pass) process.exitCode=1; };
T(1,0,'14/3',true); T(1,0,'4.6667',false); T(1,0,'28/6',true);
T(2,0,'7/4',true); T(2,1,'-3-2√3',true); T(2,1,'-3+2√3',false);
T(3,0,'505, 251, 35, 1',true); T(3,0,'1,35,251',false);
T(4,0,'(30/11, 48/11)',true); T(4,0,'(48/11, 30/11)',false);
T(5,0,'18:8',true); T(5,0,'4:9',false);
T(6,0,'a<-3, a>3',true); T(6,0,'a<-3 または a>=3',false); T(6,0,'-3<a<3',false); T(6,0,'a < -3 or 3 < a',true);
T(7,0,{mode:'cases',rows:[{cond:'m<1/2',val:'6-2m'},{cond:'m>=1/2',val:'(14m+8)/3'}]},true);
T(7,0,{mode:'cases',rows:[{cond:'m<1/2',val:'6-2m'}]},false);
T(7,0,{mode:'cases',rows:[{cond:'m<1',val:'6-2m'},{cond:'m>=1',val:'(14m+8)/3'}]},false);
T(7,0,'6-2m',false);
T(8,0,'',null);
T(9,0,'18',true); T(9,0,'18個',true);
T(10,0,'-16<k<9',true); T(10,0,'-16<=k<9',false); T(10,0,'k>-16 かつ k<9',true);
T(11,0,'4±√5',true); T(11,0,'4-√5, 4+√5',true);
T(12,0,'120',true);
T(13,0,'y=9x-16',true); T(13,0,'9x-16',true);
T(14,0,'3·2^(n-1)-n-1',true); T(14,0,'(3/2)2^n-n-1',true); T(14,0,'3·2^n-n-1',false);
T(15,0,'なし',true); T(15,0,'x>1',false);
T(16,0,{mode:'cases',rows:[{cond:'a<=0',val:'a'},{cond:'0<a<2',val:'a-a^2'},{cond:'a>=2',val:'4-3a'}]},true);
// v0.6：度の記号つきの範囲・3乗根
const one = (spec, inp, want) => { const f = parseSpec(spec).parts[0].fields[0]; const r = checkField(f, { mode: 'plain', text: inp }); if (r.ok !== want) { console.log('FAIL', spec, inp, r); process.exitCode = 1; } else passed++; };
one('(1) 30<=θ<=150', '30°<=θ<=150°', true);
one('(1) {60°, 180°}', '60, 180', true);
one('4^(1/3)-2', '2^(2/3)-2', true);
// CELESTIA v0.4：数字の添字は数字だけ（r_1v_1、B_2d）、物理量は正の値で比べる
one('(1) v_2：r_1 v_1/r_2', 'r_1v_1/r_2', true);
one('(1) R：m V/(q B B_2 d)', 'mV/(qBB_2d)', true);
one('(1) V：m sqrt(2g h/(M(M+m)))', '(m/M)√(2Mgh/(M+m))', true);
one('(1) U：-G M m/R', 'GMm/R', false);
// CELESTIA v0.6：とても小さい・大きい値も相対誤差で比べる（以前は 1e-7 の絶対誤差で、10^-11 の答えに 0 や半分の値まで正解になっていた）
one('(1) 6.24×10^-11 [J]', '6.24×10^-11', true);
one('(1) 6.24×10^-11 [J]', '6.240*10^(-11)', true);
one('(1) 6.24×10^-11 [J]', '624*10^(-13)', true);
one('(1) 6.24×10^-11 [J]', '3.12×10^-11', false);
one('(1) 6.24×10^-11 [J]', '0', false);
one('(1) 6.24×10^-11 [J]', '6.25×10^-11', false);
one('(1) 6.6×10^-34', '6.6×10^-33', false);
one('(1) 3.0×10^8 [m/s]', '3.0×10^8', true);
one('(1) 3.0×10^8 [m/s]', '3.0001×10^8', false);
one('(1) 0', '0', true);
one('(1) 0', '1-1', true);
one('(1) 0', '10^-20', false);
one('(1) E：h c/λ', 'hc/λ', true);
one('(1) E：h c/λ', '2hc/λ', false);
one('(1) F：k q^2/r^2', 'kq^2/r^2', true);
one('(1) F：k q^2/r^2', 'kq^2/(2r^2)', false);
one('(1) x：v_0 t-g t^2/2', 'v_0t-gt^2/2', true);
one('(1) x：v_0 t-g t^2/2', 'v_0t-gt^2', false);
// 全問の自動採点の答え：正解は正解、半分・0 は不正解（数値の答え）／2倍は不正解（式の答え）
{
  const fs = await import('fs'); const path = await import('path');
  const dir = path.join(path.dirname(new URL(import.meta.url).pathname), '..', 'data');
  let n = 0, bad = 0;
  for (const f of fs.readdirSync(dir).filter(f => f.startsWith('f-'))) {
    for (const [id, b] of Object.entries(JSON.parse(fs.readFileSync(path.join(dir, f), 'utf8')))) {
      if (!b.in) continue;
      let sp; try { sp = parseSpec(b.in); } catch { continue; }
      for (const pt of sp.parts) for (const fl of pt.fields) {
        if (fl.kind !== 'expr') continue;
        const src = fl.value.src, c = (t) => checkField(fl, { mode: 'plain', text: t }).ok;
        const constant = !/[a-zA-Zα-ωΑ-Ω]/.test(src.replace(/sqrt|sin|cos|tan|log|exp|π|pi/g, ''));
        const v = fl.value.fn({});
        n++;
        const zero = constant && v === 0;
        const okRight = c(src), okHalf = !zero && c(`(${src})/2`), okZero = constant && Number.isFinite(v) && v !== 0 ? c('0') : false;
        const okDouble = !zero && c(`2(${src})`);
        if (!okRight || okHalf || okZero || okDouble) { bad++; if (bad <= 5) console.log('FAIL 全問', id, src, { okRight, okHalf, okZero, okDouble }); }
      }
    }
  }
  if (bad) process.exitCode = 1; else passed++;
  console.log(`  全問の式の答え ${n}欄：正解は正解・半分／2倍／0 は不正解（だめ ${bad}）`);
}
// toTeX：TeX の命令のすぐ後に文字が続くとき（4\pik_{0} → 4\pi k_{0}）は空白を入れる（Lemma v0.28 の修正。添字・ギリシャ文字でも）
{
  const { exprTeX } = await import('../js/answer.js');
  const { createRequire } = await import('module');
  const katex = createRequire(import.meta.url)('../vendor/katex.min.js');
  for (const [src, want] of [['4πk_0', '4\\pi k_{0}'], ['2πfL', '2\\pi fL'], ['θx', '\\theta x'], ['λ_0x', '\\lambda_{0}x'], ['2πmv/(qB)', null], ['4πε_0r^2', null], ['μ_0nI', null], ['ωLt', null]]) {
    const t = exprTeX(src);
    let rendered = true; try { katex.renderToString(t, { throwOnError: true }); } catch { rendered = false; }
    if ((want && t !== want) || !rendered) { console.log('FAIL toTeX', src, t, rendered); process.exitCode = 1; } else passed++;
  }
}
if (!fieldTeX(P(7)[0]).includes('cases')) { console.log('FAIL fieldTeX'); process.exitCode = 1; }
console.log(`answer.test: ${passed} 件 OK${process.exitCode ? '（失敗あり）' : ''}`);
