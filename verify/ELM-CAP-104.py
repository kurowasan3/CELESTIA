# ELM-CAP-104: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-CAP-104.py）
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

eps0, Sa, d, V = S('eps0'), S('S'), S('d'), S('V')
C0 = eps0 * Sa / d
# 金属板：残り真空 d/2（金属内は E=0）
C1 = eps0 * Sa / (d - d / 2)
check('1', C1)
WB = (C1 * V - C0 * V) * V
check('2', WB)
WF = (C1 * V**2 / 2 - C0 * V**2 / 2) - WB
check('3', WF); assert sp.simplify(WF) .is_negative
# 誘電体（比誘電率3，厚さ d/2）：直列
C2 = 1 / ((d / 2) / (eps0 * Sa) + (d / 2) / (3 * eps0 * Sa))
Q = C0 * V
check('4', Q / C2)
U0, U2 = Q**2 / (2 * C0), Q**2 / (2 * C2)
check('5', U2 / U0)
check('6', U2 - U0)
print('OK', PID)
