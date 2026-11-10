# ELM-FLD-303: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-FLD-303.py）
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
k, q, a, Q = S('k'), S('q'), S('a'), S('Q')
A = sp.Matrix([0, a * sp.sqrt(3) / 2]); B = sp.Matrix([-a / 2, 0]); C = sp.Matrix([a / 2, 0])
G = (A + B + C) / 3
r = sp.sqrt(((A - G).T * (A - G))[0])
check('1', r, 'r')
def field(qi, P):
    v = G - P; d = sp.sqrt((v.T * v)[0]); return k * qi * v / d**3
E = field(-q, A) + field(q, B) + field(q, C)
assert sp.simplify(E[0]) == 0 and sp.simplify(E[1]) > 0   # A に向かう向き（+y）
check('2', sp.simplify(E[1]), 'E_G')
V = sum(k * qi / sp.sqrt(((G - P).T * (G - P))[0]) for qi, P in ((-q, A), (q, B), (q, C)))
check('3', sp.simplify(V), 'V_G')
check('4', sp.simplify(Q * V), 'W')
U = k * ((-q) * q + (-q) * q + q * q) / a
check('5', U, 'U')
print('OK', PID)
