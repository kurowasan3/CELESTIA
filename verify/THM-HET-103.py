# THM-HET-103 の答えを独立に導いて problem.json の入力用の答えと照合する（python3 -I verify/THM-HET-103.py）
# --- 共通：problem.json の入力用の答え（grading）を sympy の式に読み直して照合する
import json, re, sys
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
GREEK = {'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'ρ': 'rho', 'Δ': 'Delta', 'θ': 'theta', 'π': 'pi'}
TOK = re.compile(r'\s*(sqrt|\d+(?:\.\d+)?|[A-Za-zα-ωΔ](?:_(?:\{[^}]*\}|[A-Za-z0-9]+))?|\*\*|[-+*/^()])')

def to_sympy(src):
    s = src.strip()
    s = re.sub(r'\s*\[[^\]]*\]\s*$', '', s)      # 単位
    s = re.sub(r'倍$', '', s).strip()
    toks, pos = [], 0
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m: raise ValueError(f'読めない：{s[pos:]}')
        toks.append(m.group(1)); pos = m.end()
    out, names = [], {}
    def is_val(t): return t not in ('sqrt', '**', '^') and (t == ')' or re.match(r'[\dA-Za-zα-ωΔ]', t))
    def starts_val(t): return t in ('(', 'sqrt') or re.match(r'[\dA-Za-zα-ωΔ]', t)
    prev = None
    for t in toks:
        if prev is not None and is_val(prev) and starts_val(t): out.append('*')
        if t == '^': out.append('**')
        elif re.match(r'[A-Za-zα-ωΔ]', t) and t != 'sqrt':
            base = t.replace('{', '').replace('}', '')
            py = ''.join(GREEK.get(c, c) for c in base)
            names[py] = sp.Symbol(base.replace('_', '_'), positive=True)
            if py == 'pi': names[py] = sp.pi
            out.append(py)
        else: out.append(t)
        prev = t
    return sp.sympify(''.join(out), locals={**names, 'sqrt': sp.sqrt}, rational=True)

def fields(pid):
    p = json.loads((ROOT / 'content/problems/THM' / pid / 'problem.json').read_text(encoding='utf-8'))
    res = {}
    for a in p['answer']['parts']:
        g = a.get('grading')
        if not g: continue
        res[a['id']] = [to_sympy(re.sub(r'^[^：]{1,80}：', '', x.strip())) for x in g.split(' ; ')]
    return p, res

def S(name): return sp.Symbol(name, positive=True)

def check(got, want, label):
    got, want = sp.sympify(got), sp.sympify(want)
    d = sp.simplify(got - want)
    ok = d == 0
    if not ok and d.free_symbols:   # 平方根の整理の違いなど：正の乱数を何組か代入して確かめる
        import random
        rnd = random.Random(0)
        ok = all(abs(complex(d.subs({s_: rnd.uniform(0.5, 3) for s_ in d.free_symbols}).evalf())) < 1e-9 for _ in range(20))
    print(('OK  ' if ok else 'NG  ') + f'{label}: problem.json = {got}   derived = {sp.simplify(want)}')
    if not ok: sys.exit(f'食い違い：{label}')

# --- この問題
p, g = fields('THM-HET-103')
cw, ci, Lf = sp.Rational(42, 10), sp.Rational(21, 10), 336
mw, Tw = 200, 30

def mix(mi, Ti):
    """氷 mi [g]（温度 Ti）を入れたときの (最終温度, とけた質量)。場合分けを機械的に行う"""
    warm = mi * ci * (0 - Ti)                 # 氷を 0℃ にする熱
    avail = mw * cw * Tw                     # 水が 0℃ まで冷えて出せる熱
    assert avail > warm                      # 水が凍り始めることはない（この問題の範囲）
    if avail >= warm + mi * Lf:              # すべてとける
        t = sp.Symbol('t')
        T = sp.solve(sp.Eq(mw * cw * (Tw - t), warm + mi * Lf + mi * cw * t), t)[0]
        assert T >= 0
        return T, mi
    x = (avail - warm) / Lf
    assert 0 <= x <= mi
    return sp.Integer(0), x

check(g['1'][0], 50 * ci * 20, '(1) 氷を 0℃ にする熱 [J]')
T1, x1 = mix(50, -20)
assert x1 == 50, 'すべてとけるはず'
check(g['3'][0], T1, '(3) 最終温度 [℃]')
T2, x2 = mix(100, -20)
check(g['4'][0], T2, '(4) 最終温度 [℃]')
check(g['4'][1], x2, '(4) とけた質量 [g]')
print('THM-HET-103: すべて一致')
