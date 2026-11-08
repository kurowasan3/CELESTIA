# ELM-DCC-104: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-DCC-104.py）
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

V = sp.Symbol('V', positive=True)
Ib = lambda v: sp.Rational(1, 10) * sp.sqrt(v)
R = 10
v1 = sp.solve(sp.Eq(6, V + R * Ib(V)), V); assert len(v1) == 1; v1 = v1[0]
check('1', v1, 'V'); check('1', Ib(v1), 'I')
v2 = sp.solve(sp.Eq(6, 2 * V + R * Ib(V)), V); assert len(v2) == 1; v2 = v2[0]
check('2', v2, 'V'); check('2', Ib(v2), 'I')
check('3', v1 / Ib(v1), 'R_1'); check('3', v2 / Ib(v2), 'R_2')
assert v1 * Ib(v1) > v2 * Ib(v2)   # 電力が大きい方が高温・高抵抗
v4 = sp.solve(sp.Eq(V * Ib(V), sp.Rational(1, 10)), V)[0]
check('4', v4 + R * Ib(v4), 'E')
# 図の曲線が式どおり（V=6 で約 0.245 A）
assert abs(float(Ib(6)) - 0.245) < 1e-3
print('OK', PID)
