# THM-FST-205 の答えを独立に導いて problem.json の入力用の答えと照合する（python3 -I verify/THM-FST-205.py）
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
p, g = fields('THM-FST-205')
R = sp.Rational
P0, V0, T0, n, Rg = S('P_0'), S('V_0'), S('T_0'), S('n'), S('R')
gam = (R(5, 2)) / R(3, 2)
V1 = V0 / 8
P1 = P0 * (V0 / V1) ** gam
T1 = T0 * (P1 * V1) / (P0 * V0)
check(g['1'][0], P1, '(1) P1')
check(g['2'][0], T1, '(2) T1')
# 外部がした仕事 = -∫P dV（断熱曲線に沿って数値でなく記号で積分）
V = sp.Symbol('V', positive=True)
Wad = -sp.integrate(P0 * (V0 / V) ** gam, (V, V0, V1))
check(g['3'][0], sp.simplify(Wad), '(3) 断熱圧縮の仕事（∫PdV）')
Wiso = -sp.integrate(P0 * V0 / V, (V, V0, V1))
check(g['4'][0], sp.simplify(Wiso), '(4) 等温圧縮の仕事（∫PdV）')
assert float((Wad - Wiso).subs({P0: 1, V0: 1})) > 0
check_sig(g['5'][0], 300 * 16 ** 0.4, '(5) ディーゼル圧縮後の温度 [K]')
check_sig(g['5'][0], 300 * 3.03, '(5) （16^0.4=3.03 を使った値）')
print('THM-FST-205: すべて一致')
