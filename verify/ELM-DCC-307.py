# ELM-DCC-307: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-DCC-307.py）
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
x, E, R0, RL = S('x'), S('E'), S('R_0'), S('R_L')
check('1', x * E, 'V')
Rpb = x * R0 * RL / (x * R0 + RL)
V = sp.simplify(E * Rpb / (Rpb + (1 - x) * R0)); check('2', V, 'V')
assert sp.limit(V, RL, sp.oo) == x * E
check('3', V.subs({x: R_(1, 2), RL: R0}), 'V')
ratio = sp.simplify(V / (x * E))
# x(1-x) の最大は x=1/2 → 条件 ratio(1/2) = 0.99 となる R_L
RLmin = sp.solve(sp.Eq(ratio.subs(x, R_(1, 2)), R_(99, 100)), RL)[0]
check('4', RLmin, 'R_L')
import numpy as np
xs = np.linspace(0.001, 0.999, 999); f = sp.lambdify((x, R0, RL), ratio)
assert f(xs, 1.0, float(RLmin.subs(R0, 1))).min() >= 0.99 - 1e-12
Itot = E / ((1 - x) * R0 + Rpb)
eta = sp.simplify((V**2 / RL) / (E * Itot)).subs({x: R_(1, 2), RL: R0})
check('5', sp.simplify(eta), 'η')
print('OK', PID)
