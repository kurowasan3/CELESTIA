# ELM-FLD-103: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-FLD-103.py）
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

X = sp.Symbol('x')
V = sp.Piecewise((300 * X, X <= sp.Rational(1, 10)), (30, X <= sp.Rational(2, 10)), (30 * (sp.Rational(4, 10) - X) / sp.Rational(2, 10), X <= sp.Rational(4, 10)), (0, True))
Ex = lambda x0: -sp.diff(V, X).subs(X, x0)
assert Ex(sp.Rational(5, 100)) == -300   # 負の向き
check('1', abs(Ex(sp.Rational(5, 100))))
assert Ex(sp.Rational(3, 10)) == 150 and Ex(sp.Rational(15, 100)) == 0
check('2', Ex(sp.Rational(3, 10)))
e = sp.Rational(16, 10) * sp.Rational(1, 10**19)
K = lambda x0: (-e) * V.subs(X, 0) - (-e) * V.subs(X, x0)   # 運動エネルギー（エネルギー保存）
check('3', K(sp.Rational(15, 100)))
# 最大到達点：x>0 で K=0 となる最初の点（数値で走査）
xs = [sp.Rational(i, 1000) for i in range(1, 500)]
xmax = next(x0 for x0 in xs if K(x0) <= 0)
check('4', xmax)
print('OK', PID)
