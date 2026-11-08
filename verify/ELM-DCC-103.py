# ELM-DCC-103: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-DCC-103.py）
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

I1, I2, I3 = sp.symbols('I1 I2 I3')
E1, E2, R1, R2, R3 = 6, 4, 2, 1, 2
sol = sp.solve([I1 + I2 - I3, E1 - R1 * I1 - R3 * I3, E2 - R2 * I2 - R3 * I3], [I1, I2, I3])
check('1', sol[I1]); check('2', sol[I2]); check('3', sol[I3])
check('4', R3 * sol[I3])
check('5', R3 * sol[I3]**2)
assert E1 * sol[I1] + E2 * sol[I2] == R1 * sol[I1]**2 + R2 * sol[I2]**2 + R3 * sol[I3]**2
fig = [f for f in P['figures'] if f['id'] == 'a-currents'][0]['physics']['junctions'][0]['values']
assert abs(fig['current:I1'] - sol[I1]) < 1e-9 and abs(fig['current:I2'] - sol[I2]) < 1e-9 and abs(fig['current:I3'] - sol[I3]) < 1e-9
print('OK', PID)
