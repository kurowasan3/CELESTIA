# ELM-DCC-308: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-DCC-308.py）
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
E, r, R = 12.0, 1.0, 2.0
I = E / (R + r); P = R * I * I
check('1', sp.Float(I), 'I'); check('1', sp.Float(P), 'P')
slope = (60 - 20) / 2000         # 図の接線
Ch = P / slope; check('2', sp.Float(Ch), 'C_h')
k = P / (60 - 20); check('3', sp.Float(k), 'k')
check('4', sp.Float((P - k * 20) / Ch))
# 数値で熱の収支を積分し，図の曲線（20 + 40(1 - e^{-t/2000})）と一致することを確かめる
T, dt = 20.0, 0.1
for i in range(1, 60001):
    T += (P - k * (T - 20)) / Ch * dt
    if i % 10000 == 0:
        tt = i * dt; assert abs(T - (20 + 40 * (1 - np.exp(-tt / 2000)))) < 1e-3
Rs = np.linspace(0.1, 5, 49001); Ps = E**2 * Rs / (Rs + r)**2
Rb = Rs[np.argmax(Ps)]
check('5', sp.Float(round(Rb, 6)), 'R'); check('5', sp.Float(20 + Ps.max() / k), 'T')
Chs, ks, Ps_ = S('C_h'), S('k'), S('P')
t1 = sp.solve(sp.Eq(Ps_ / Chs * sp.Symbol('t'), Ps_ / ks), sp.Symbol('t'))[0]
check('6', t1, 't_1')
assert abs(float(t1.subs({Chs: Ch, ks: k})) - 2000) < 1e-9
print('OK', PID)
