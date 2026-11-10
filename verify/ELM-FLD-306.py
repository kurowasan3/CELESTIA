# ELM-FLD-306: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-FLD-306.py）
# --- 共通：problem.json の grading を読み，sympy の式にして独立に導いた答えと照合する
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor
PID = Path(__file__).stem
P = json.loads((Path(__file__).resolve().parent.parent / 'content/problems/ELM' / PID / 'problem.json').read_text(encoding='utf-8'))
GR = {x['id']: x.get('grading') for x in P['answer']['parts']}
TR = standard_transformations + (implicit_multiplication_application, convert_xor)
GREEK = {'π': 'pi', 'ε_0': 'eps0', 'ε_r': 'epsr', 'θ_0': 'theta0', 'θ': 'theta'}
def to_expr(s):
    s = re.sub(r'\s*\[[^\[\]]+\]\s*$', '', s.strip())
    s = re.sub(r'倍$|回$', '', s).strip()
    for g, n in GREEK.items():
        s = s.replace(g, n)
    s = re.sub(r'_\{(\w+)\}', r'_\1', s)
    names = set(re.findall(r'[A-Za-z][A-Za-z0-9_]*', s)) - {'sqrt', 'pi'}
    loc = {n: sp.Symbol(n, positive=True) for n in names}
    return parse_expr(s, local_dict=loc, transformations=TR)
def fields(pid):
    out = {}
    for i, f in enumerate([x.strip() for x in GR[pid].split(' ; ')]):
        lab, _, body = f.rpartition('：')
        out[lab or i] = to_expr(body)
    return out
def same(a, b, rel=1e-9):
    a, b = sp.sympify(a), sp.sympify(b)
    syms = sorted(a.free_symbols | b.free_symbols, key=str)
    for _ in range(5):
        sub = {x: sp.Rational(random.randint(11, 97), random.randint(5, 13)) for x in syms}
        va, vb = sp.N(a.subs(sub), 30), sp.N(b.subs(sub), 30)
        if abs(va - vb) > rel * max(1, abs(vb)):
            return False
    return True
def check(pid, expected, label=None, sig=None):
    """sig を与えると，導いた値を有効数字 sig 桁に丸めたものと grading が一致するか確かめる"""
    fs = fields(pid)
    got = fs[label] if label is not None else list(fs.values())[0]
    if sig is not None:
        x = float(sp.N(expected)); expected = sp.Float(f'{x:.{sig - 1}e}')
    assert same(got, expected), f'({pid}{":" + str(label) if label is not None else ""}) grading {got} != derived {sp.N(expected)}'
    print(f'  ({pid}{":" + str(label) if label is not None else ""}) {got}  ✓')
S = lambda n: sp.Symbol(n, positive=True)
R_ = sp.Rational
# ---
import numpy as np
xs = np.linspace(0, 0.45, 450001)
E = np.piecewise(xs, [xs <= 0.10, (xs > 0.10) & (xs <= 0.30), (xs > 0.30) & (xs <= 0.40), xs > 0.40],
                 [2000, lambda x: 2000 - 20000 * (x - 0.10), lambda x: -2000 + 20000 * (x - 0.30), 0])
V = -np.concatenate([[0], np.cumsum((E[1:] + E[:-1]) / 2 * np.diff(xs))])
at = lambda x: V[int(round(x / 0.45 * 450000))]
check('1', sp.Float(round(at(0.10), 3)))
i = np.argmin(V[xs <= 0.40]); check('2', sp.Float(round(xs[i], 4)), 'x'); check('2', sp.Float(round(V[i], 3)), 'V')
check('3', sp.Float(round(at(0.40), 3)))
e = 1.6e-19
# 陽子：x=0.30 から，同じ電位の点まで往復
sel = (V <= at(0.30) + 1e-6) & (xs <= 0.40)
assert abs(xs[sel].min() - 0.10) < 1e-3 and abs(xs[sel].max() - 0.30) < 1e-3
check('4', sp.Float(e * (at(0.30) - V.min())), 'K')
# 電子：位置エネルギー -eV の最大（0〜0.40）
Umax = -e * V[xs <= 0.40].min()
check('5', sp.Float(Umax), 'K')
check('6', sp.Float(Umax - (-e * at(0.40))), 'K')
print('OK', PID)
