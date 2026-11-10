# THM-FST-207 の答えを独立に導いて problem.json の入力用の答えと照合する（python3 -I verify/THM-FST-207.py）
# --- 共通：problem.json の入力用の答え（grading）を sympy の式に読み直して照合する
import json, re, sys
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
GREEK = {'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'ρ': 'rho', 'Δ': 'Delta', 'θ': 'theta', 'π': 'pi', 'μ': 'mu'}
FUNCS = {'sqrt': sp.sqrt, 'log': sp.log, 'cos': sp.cos, 'sin': sp.sin}
TOK = re.compile(r'\s*(sqrt|log|cos|sin|\d+(?:\.\d+)?|[A-Za-zα-ωΔ](?:_(?:\{[^}]*\}|[A-Za-z0-9]+))?|\*\*|[-+*/^()])')


def to_sympy(src):
    s = src.strip()
    s = re.sub(r'\s*\[[^\]]*\]\s*$', '', s)      # 単位
    s = re.sub(r'倍$', '', s).strip()
    toks, pos = [], 0
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m:
            raise ValueError(f'読めない：{s[pos:]}')
        toks.append(m.group(1)); pos = m.end()
    out, names = [], {}
    def is_val(t): return t not in FUNCS and t not in ('**', '^') and (t == ')' or re.match(r'[\dA-Za-zα-ωΔ]', t))
    def starts_val(t): return t == '(' or t in FUNCS or re.match(r'[\dA-Za-zα-ωΔ]', t)
    prev = None
    for t in toks:
        if prev is not None and is_val(prev) and starts_val(t):
            out.append('*')
        if t == '^':
            out.append('**')
        elif t in FUNCS:
            out.append(t)
        elif re.match(r'[A-Za-zα-ωΔ]', t):
            base = t.replace('{', '').replace('}', '')
            py = ''.join(GREEK.get(c, c) for c in base)
            names[py] = sp.pi if py == 'pi' else sp.Symbol(base, positive=True)
            out.append(py)
        else:
            out.append(t)
        prev = t
    return sp.sympify(''.join(out), locals={**names, **FUNCS}, rational=True)


def fields(pid):
    p = json.loads((ROOT / 'content/problems/THM' / pid / 'problem.json').read_text(encoding='utf-8'))
    res = {}
    for a in p['answer']['parts']:
        g = a.get('grading')
        if not g:
            continue
        res[a['id']] = [to_sympy(re.sub(r'^[^：]{1,80}：', '', x.strip())) for x in g.split(' ; ')]
    return p, res


def S(name):
    return sp.Symbol(name, positive=True)


def check(got, want, label, rel=None):
    """rel を与えると相対誤差 rel 以内（有効数字で丸めた値の照合）"""
    got, want = sp.sympify(got), sp.sympify(want)
    if rel is not None:
        g, w = float(got), float(want)
        ok = abs(g - w) <= rel * abs(w)
    else:
        d = sp.simplify(got - want)
        ok = d == 0
        if not ok and d.free_symbols:
            import random
            rnd = random.Random(0)
            ok = all(abs(complex(d.subs({s_: rnd.uniform(0.5, 3) for s_ in d.free_symbols}).evalf())) < 1e-9 for _ in range(20))
        elif not ok:
            ok = abs(complex(d.evalf())) < 1e-12
    print(('OK  ' if ok else 'NG  ') + f'{label}: problem.json = {got}   derived = {sp.N(want, 8) if rel is not None else sp.simplify(want)}')
    if not ok:
        sys.exit(f'食い違い：{label}')


def sig(x, n=2):
    """有効数字 n 桁に丸めた値（問題文で有効数字を指定した答えの照合用）"""
    x = float(x)
    from math import floor, log10
    k = floor(log10(abs(x)))
    return round(x, -k + n - 1)


def check_sig(got, want, label, n=2):
    """problem.json の値が，導いた値を有効数字 n 桁に丸めたものと一致するか"""
    r = sig(want, n)
    ok = abs(float(got) - r) <= 1e-9 * abs(r)
    print(('OK  ' if ok else 'NG  ') + f'{label}: problem.json = {got}   derived = {float(want):.6g} → 有効数字{n}桁 {r:g}')
    if not ok:
        sys.exit(f'食い違い：{label}')

# --- この問題
p, g = fields('THM-FST-207')
R = sp.Rational
n, Rg, T0, P0 = S('n'), S('R'), S('T_0'), S('P_0')
gam = R(5, 3)
# 状態を V と P から独立に求める：A(P0, VA)，断熱 PV^γ 一定
VA = n * Rg * T0 / P0
PB = 32 * P0
VB = VA * (P0 / PB) ** (1 / gam)
TB = PB * VB / (n * Rg)
TC = 16 * T0
VC = n * Rg * TC / PB
VD = VC * (PB / P0) ** (1 / gam)
TD = P0 * VD / (n * Rg)
check(g['1'][0], TB, '(1) T_B')
check(g['2'][0], TD, '(2) T_D')
Qin = R(5, 2) * n * Rg * (TC - TB)
Qout = R(5, 2) * n * Rg * (TD - T0)
check(g['3'][0], Qin, '(3) Q_in')
check(g['3'][1], Qout, '(3) Q_out')
# 仕事を面積（∫PdV）で独立に：2つの定圧 + 2つの断熱
V = sp.Symbol('V', positive=True)
W = PB * (VC - VB) + sp.integrate(PB * (VC / V) ** gam, (V, VC, VD)) + P0 * (VA - VD) + sp.integrate(P0 * (VA / V) ** gam, (V, VA, VB))
check(g['4'][0], sp.simplify(W), '(4) W（∫PdV）')
check(g['4'][1], sp.simplify(W / Qin), '(4) e')
r = S('r')
TBr = T0 * r ** ((gam - 1) / gam)
TDr = TC / r ** ((gam - 1) / gam)
e_r = 1 - (TDr - T0) / (TC - TBr)
check(g['5'][0], sp.simplify(e_r), '(5) e(r)')
x = sp.Symbol('x', positive=True)
Wx = R(5, 2) * n * Rg * ((TC - x) - (T0 * TC / x - T0))
xs = [s_ for s_ in sp.solve(sp.diff(Wx, x), x) if s_.is_positive]
check(g['6'][0], xs[0], '(6) W 最大の T_B')
assert sp.simplify(sp.diff(Wx, x, 2).subs(x, xs[0])) < 0 if xs[0].free_symbols == set() else True
print('    (7) T_B→16T0 で e →', sp.simplify(1 - T0 / TC), '，W →', sp.simplify(Wx.subs(x, TC)))
assert sp.simplify(Wx.subs(x, TC)) == 0
print('THM-FST-207: すべて一致')
