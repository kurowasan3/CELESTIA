# ELM-FLD-307: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-FLD-307.py）
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
m, g, q, E, l, th0 = S('m'), S('g'), S('q'), S('E'), S('l'), S('theta0')
check('1', q * E / (m * g), 'tanθ_0')
S0 = sp.sqrt((m * g)**2 + (q * E)**2); check('2', S0, 'S_0')
g1 = S0 / m
check('3', 2 * sp.pi * sp.sqrt(l / g1), 'T')
check('4', 2 * th0, 'θ_max')
# 数値で：真下から放した運動を積分し，最大角・最大速さ・そのときの張力を確かめる
vals = dict(m=0.2, g=9.8, q=1e-3, E=1500.0, l=0.8)
mm, gg, qq, EE, ll = (vals[k] for k in 'm g q E l'.split())
t0 = np.arctan(qq * EE / (mm * gg))
def acc(t): return (qq * EE * np.cos(t) - mm * gg * np.sin(t)) / (mm * ll)
th, w, dt, thmax, wmax, thw = 0.0, 0.0, 1e-5, 0.0, 0.0, 0.0
for _ in range(400000):
    a1 = acc(th); w += a1 * dt / 2; th += w * dt; w += acc(th) * dt / 2   # 速度ベルレ
    thmax = max(thmax, th)
    if abs(w) > wmax: wmax, thw = abs(w), th
    if w < 0 and th < 0.01: break
assert abs(thmax - 2 * t0) < 1e-3, (thmax, 2 * t0)
assert abs(thw - t0) < 1e-3
sub = {m: mm, g: gg, q: qq, E: EE, l: ll}
vmax = fields('5')['v_max'].subs(sub)
assert abs(float(vmax) - wmax * ll) < 1e-3 * wmax * ll
# 張力：糸方向の運動方程式 S - (mg cosθ + qE sinθ) = m v^2 / l（θ = θ0）
S1 = mm * gg * np.cos(t0) + qq * EE * np.sin(t0) + mm * (wmax * ll)**2 / ll
assert abs(float(fields('6')['S_1'].subs(sub)) - S1) < 1e-3 * S1
check('5', sp.sqrt(2 * l * (g1 - g)), 'v_max'); check('6', 3 * S0 - 2 * m * g, 'S_1')
# 張力が負にならない：S(φ) = m g1 (3cosφ - 2cosθ0) ≥ m g1 cosθ0 > 0（|φ|≤θ0）
print('OK', PID)
