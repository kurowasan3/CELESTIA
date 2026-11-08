# ELM-FLD-104: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-FLD-104.py）
# --- 共通：problem.json の grading を読み，sympy の式にして独立に導いた答えと照合する
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor
PID = Path(__file__).stem
P = json.loads((Path(__file__).resolve().parent.parent / 'content/problems/ELM' / PID / 'problem.json').read_text(encoding='utf-8'))
GR = {x['id']: x.get('grading') for x in P['answer']['parts']}
TR = standard_transformations + (implicit_multiplication_application, convert_xor)
def to_expr(s):
    s = re.sub(r'\s*\[[^\[\]]+\]\s*$', '', s.strip())
    s = re.sub(r'倍$|回$', '', s).strip()
    s = s.replace('π', 'pi').replace('ε_0', 'eps0')
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
def same(a, b):
    a, b = sp.sympify(a), sp.sympify(b)
    syms = sorted(a.free_symbols | b.free_symbols, key=str)
    for _ in range(5):
        sub = {x: sp.Rational(random.randint(11, 97), random.randint(5, 13)) for x in syms}
        va, vb = sp.N(a.subs(sub), 30), sp.N(b.subs(sub), 30)
        if abs(va - vb) > 1e-9 * max(1, abs(vb)):
            return False
    return True
def check(pid, expected, label=None):
    fs = fields(pid)
    got = fs[label] if label is not None else list(fs.values())[0]
    assert same(got, expected), f'({pid}) grading {got} != derived {sp.nsimplify(expected)}'
    print(f'  ({pid}{":" + str(label) if label is not None else ""}) {got}  ✓')
S = lambda n: sp.Symbol(n, positive=True)
# ---

k, Q, a, q, m, y = S('k'), S('Q'), S('a'), S('q'), S('m'), S('y')
V = 2 * k * Q / sp.sqrt(a**2 + y**2)
E = sp.simplify(-sp.diff(V, y))
check('1', E, 'E')
crit = [c for c in sp.solve(sp.diff(E, y), y) if c.is_positive]
assert len(crit) == 1
assert sp.simplify(sp.diff(E, y, 2).subs(y, crit[0])) < 0
check('2', crit[0], 'y')
vinf = sp.sqrt(2 * q * V.subs(y, a) / m)
check('3', vinf, 'v_∞')
vO = sp.sqrt(2 * (-q * V.subs(y, a) + q * V.subs(y, 0)) / m)
check('4', vO, 'v_O')
Ky = sp.limit(q * E / y, y, 0)   # 力 -qE ≈ -Ky*y
T = 2 * sp.pi * sp.sqrt(m / Ky)
check('5', T, 'T')
# 数値で単振動の周期を確かめる（微小振幅）
import numpy as np
kv, Qv, av, qv, mv = 1.0, 1.0, 1.0, 1.0, 1.0
f = lambda yy: -qv * 2 * kv * Qv * yy / (av**2 + yy**2)**1.5 / mv
yy, vv, t, dt, cross = 1e-3, 0.0, 0.0, 1e-4, []
while len(cross) < 3:
    a1 = f(yy); vv += a1 * dt / 2; yn = yy + vv * dt; vv += f(yn) * dt / 2
    if yy > 0 >= yn: cross.append(t)
    yy = yn; t += dt
assert abs((cross[2] - cross[1]) - float(T.subs({k: 1, Q: 1, a: 1, q: 1, m: 1}))) < 1e-3
print('OK', PID)
