# ELM-CAP-103: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-CAP-103.py）
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

import numpy as np
E, I0, tau = sp.Rational(6), sp.Rational(3, 1000), sp.Rational(2, 10)
R = E / I0; C = tau / R
check('1', R / 1000); check('2', C * 10**6)
# 数値積分でグラフの面積 = Q を確かめる
t = np.linspace(0, 20 * float(tau), 200001); I = float(I0) * np.exp(-t / float(tau))
area = np.trapz(I, t) if hasattr(np, 'trapz') else np.trapezoid(I, t)
assert abs(area - float(C * E)) < 1e-8
check('3', C * E, 'Q')
check('4', (E - 2) / R * 1000)
heat = np.trapz(I**2 * float(R), t) if hasattr(np, 'trapz') else np.trapezoid(I**2 * float(R), t)
assert abs(heat - float(C * E**2 / 2)) < 1e-8
check('5', C * E**2 / 2)
# 接線の t 切片 = τ
assert sp.simplify(I0 / (I0 / tau) - tau) == 0
print('OK', PID)
