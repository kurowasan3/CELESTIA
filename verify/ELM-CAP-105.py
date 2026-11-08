# ELM-CAP-105: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-CAP-105.py）
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

C1, C2, E, n = S('C_1'), S('C_2'), S('E'), S('n')
V1 = C1 * E / (C1 + C2)
check('1', V1, 'V_1')
check('2', C1 * E**2 / 2 - (C1 + C2) * V1**2 / 2, 'J_1')
# 漸化式を数値で回して一般項と比べる
vals = {C1: sp.Rational(3), C2: sp.Rational(5), E: sp.Rational(7)}
Vn_formula = fields('3')['V_n']
v = 0
heat = vals[C1] * vals[E]**2 / 2   # 最初の C1 の充電で出る熱（電池の仕事 C1E^2 - 蓄積 C1E^2/2）
for k in range(1, 9):
    c1, c2, e = vals[C1], vals[C2], vals[E]
    vnew = (c1 * e + c2 * v) / (c1 + c2)
    heat += (c1 * e**2 / 2 + c2 * v**2 / 2) - (c1 + c2) * vnew**2 / 2          # b 側にしたときの熱
    heat += c1 * (e - vnew) * e - (c1 * e**2 / 2 - c1 * vnew**2 / 2)              # a 側で再充電する熱
    v = vnew
    assert sp.simplify(Vn_formula.subs({**vals, n: k}) - v) == 0
print('  (3) 漸化式 n=1..8 と一致 ✓')
# C1=C2
r = sp.Rational(1, 2)
nmin = next(k for k in range(1, 50) if 1 - r**k >= sp.Rational(9, 10))
check('4', nmin)
# 熱の合計：n→∞ の極限（数値で多数回）
v = 0; c1, c2, e = vals[C1], vals[C2], vals[E]; heat = c1 * e**2 / 2
for k in range(400):
    vnew = (c1 * e + c2 * v) / (c1 + c2)
    heat += (c1 * e**2 / 2 + c2 * v**2 / 2) - (c1 + c2) * vnew**2 / 2 + c1 * (e - vnew) * e - (c1 * e**2 / 2 - c1 * vnew**2 / 2)
    v = vnew
Jf = fields('5')['J']
assert abs(float(heat) - float(Jf.subs(vals))) < 1e-9
print('  (5) 数値の和と一致 ✓')
print('OK', PID)
