# ELM-FLD-305: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-FLD-305.py）
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
e, m, V0, V, d, l, L = S('e'), S('m'), S('V_0'), S('V'), S('d'), S('l'), S('L')
v0 = sp.sqrt(2 * e * V0 / m); a = e * V / (m * d); t = l / v0
y1 = sp.simplify(a * t**2 / 2); tan = sp.simplify(a * t / v0)
check('3', y1, 'y_1'); check('4', tan, 'tanθ')
Y = sp.simplify(y1 + L * tan)
check('5', Y, 'Y')
assert not (Y.free_symbols & {e, m})       # 比電荷によらない
# 数値積分で軌道を確かめる
import numpy as np
vals = {e: 1.6e-19, m: 9.1e-31, V0: 1000, V: 50, d: 0.01, l: 0.04, L: 0.2}
vx = float(v0.subs(vals)); ay = float(a.subs(vals)); N = 100000; dt = vals[l] / vx / N; y = vy = 0.0
for _ in range(N):
    y += vy * dt + ay * dt * dt / 2; vy += ay * dt
assert abs(y - float(y1.subs(vals))) < 1e-3 * float(y1.subs(vals)) + 1e-9
# 条件 y1 < d/2
assert sp.simplify(sp.solve(sp.Eq(y1, d / 2), V)[0] - 2 * d**2 * V0 / l**2) == 0
print('OK', PID)
