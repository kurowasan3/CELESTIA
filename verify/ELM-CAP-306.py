# ELM-CAP-306: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-CAP-306.py）
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
E, R1, R2, R3, C = 12.0, 20.0, 30.0, 60.0, 50e-6
# 回路をそのまま時間発展させて確かめる（節点 P の電位 x，C の電圧 vc）
def node(vc, closed):
    # P での電流保存：(E - x)/R1 * closed = x/R2 + (x - vc)/R3
    if closed: return (E / R1 + vc / R3) / (1 / R1 + 1 / R2 + 1 / R3)
    return (vc / R3) / (1 / R2 + 1 / R3)
vc = 0.0
x = node(vc, True); I0 = (E - x) / R1; I3 = (x - vc) / R3
check('1', sp.Float(I0), 'I_0'); check('2', sp.Float(I3), 'I_3')
dt = 1e-6
for _ in range(200000):
    x = node(vc, True); vc += (x - vc) / R3 / C * dt
x = node(vc, True); I1 = (E - x) / R1
check('3', sp.Float(round(I1, 9)), 'I_1'); check('4', sp.Float(round(C * vc, 12)), 'Q')
# S を開く
x = node(vc, False); I2 = x / R2
assert I2 > 0   # P→Q（開く前と同じ向き）
check('5', sp.Float(round(I2, 9)), 'I_2')
H3 = 0.0
for _ in range(200000):
    x = node(vc, False); i = (vc - x) / R3; H3 += i * i * R3 * dt; vc -= i / C * dt
assert abs(H3 - 8.64e-4) < 2e-6, H3
check('6', R_(60, 90) * C * R_(72, 10)**2 / 2, 'H')
print('OK', PID)
