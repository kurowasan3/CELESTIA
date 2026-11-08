# ELM-DCC-105: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-DCC-105.py）
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

R1, R2, R3, R4, R5, E = [S(n) for n in ('R_1', 'R_2', 'R_3', 'R_4', 'R_5', 'E')]
x, y = sp.symbols('x y')
def solve(vals):
    r1, r2, r3, r4, r5, e = [vals[s] for s in (R1, R2, R3, R4, R5, E)]
    s = sp.solve([(e - x) / r1 - x / r2 - (x - y) / r5, (e - y) / r3 + (x - y) / r5 - y / r4], [x, y])
    return s[x], s[y], (e - s[x]) / r1 + (e - s[y]) / r3
# (1)：つり合い条件で R5 の電流が 0
R4bal = fields('1')['R_4']
for _ in range(3):
    import random as _r
    vals = {R1: _r.randint(1, 9), R2: _r.randint(1, 9), R3: _r.randint(1, 9), R5: _r.randint(1, 9), E: 12}
    vals[R4] = R4bal.subs(vals)
    xs, ys, _ = solve(vals); assert sp.simplify(xs - ys) == 0
print('  (1) つり合い条件で I5=0 ✓')
num = {R1: 2, R2: 4, R3: 4, R4: 2, R5: 2, E: 12}
Vc, Vd = sp.Rational(12 * 4, 6), sp.Rational(12 * 2, 6)
check('2', Vc - Vd)
xs, ys, Ibat = solve(num)
assert xs > ys   # C→D
check('3', (xs - ys) / 2)
check('4', Ibat)
check('5', 12 / Ibat)
# テブナンの別解とも一致
assert sp.Rational(4) / (sp.Rational(8, 3) + 2) == (xs - ys) / 2
print('OK', PID)
