# ELM-FLD-105: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-FLD-105.py）
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

k, Q, R, r = S('k'), S('Q'), S('R'), S('r')
# 球殻の電荷 q のつくる電位：r<R で kq/R，r>R で kq/r
def Vshell(qs, rr, inside): return k * qs / R if inside else k * qs / rr
# 状態1：内面 -Q，外面 +Q（ガウス：内側の閉曲面の電気量 0、全体 0）
qin, qout = -Q, Q - (-Q) - Q  # 外面 = 全電荷0 - 内面
qout = 0 - qin
f = fields('1'); assert same(f['内側'], qin) and same(f['外側'], qout); print('  (1) ✓')
V1 = k * Q / R + k * qin / R + k * qout / R
check('2', V1)
check('3', k * Q / (R / 2) + Vshell(qin, 0, True) + Vshell(qout, 0, True))
# 接地：球殻電位 0 → 外面電荷 x：kQ/R - kQ/R + kx/R = 0
x = sp.Symbol('x'); qout2 = sp.solve(sp.Eq(k * Q / R + k * qin / R + k * x / R, 0), x)[0]
check('4', qout - qout2)     # 大地へ移った正電荷
check('5', k * Q / (R / 2) + k * qin / R + k * qout2 / R)
check('6', k * (qin + qout2) / R)
print('OK', PID)
