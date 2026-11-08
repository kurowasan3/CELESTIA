# --- 共通：problem.json の入力用の答え（grading）を sympy の式に直し，独立に導いた式と数値で照合する
import json, re, random, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor
ROOT = pathlib.Path(__file__).resolve().parent.parent
TR = standard_transformations + (implicit_multiplication_application, convert_xor)


def grading(pid, part):
    f = next(ROOT.glob(f'content/problems/*/{pid}/problem.json'))
    p = json.loads(f.read_text(encoding='utf-8'))
    g = next(x for x in p['answer']['parts'] if x['id'] == part)['grading']
    out = []
    for fld in re.split(r'\s;\s', g):
        if '：' in fld:
            fld = fld.split('：', 1)[1]
        fld = re.sub(r'\[[^\]]*\]\s*$', '', fld).strip()
        fld = re.sub(r'倍$', '', fld).strip()
        s = fld.replace('π', 'pi').replace('θ', 'theta').replace('μ', 'mu').replace('ω', 'omega')
        s = re.sub(r'_\{([A-Za-z0-9]+)\}', r'_\1', s)
        names = set(re.findall(r'[A-Za-z_][A-Za-z_0-9]*', s)) - {'sqrt', 'sin', 'cos', 'tan', 'pi'}
        loc = {n: sp.Symbol(n, positive=True) for n in names}
        out.append(parse_expr(s, local_dict=loc, transformations=TR))
    return out


def same(a, b, trials=6):
    a, b = sp.sympify(a), sp.sympify(b)
    allsyms = a.free_symbols | b.free_symbols
    ok = 0
    for _ in range(trials * 20):
        if ok >= trials:
            break
        vals = {}
        for n in {s.name for s in allsyms}:
            if n == 'e':
                vals[n] = random.uniform(0.1, 0.9)
            elif n == 'theta':
                vals[n] = random.uniform(0.2, 1.3)
            elif n in ('mu', 'mu_0'):
                vals[n] = random.uniform(0.05, 0.3)
            else:
                vals[n] = random.uniform(1.0, 3.0)
        sub = {s: vals[s.name] for s in allsyms}
        x, y = complex(a.subs(sub).evalf()), complex(b.subs(sub).evalf())
        if abs(x.imag) > 1e-12 or abs(y.imag) > 1e-12:
            continue      # 物理的でない値の組（平方根の中が負）は使わない
        ok += 1
        if abs(x - y) > 1e-9 * max(1, abs(x), abs(y)):
            return False
    return ok >= trials


def check(pid, part, *expected):
    got = grading(pid, part)
    assert len(got) == len(expected), (pid, part, got)
    for g, e in zip(got, expected):
        assert same(g, e), f'{pid} ({part}): grading {g} ≠ 導出 {sp.simplify(e)}'
    print(f'{pid} ({part}) OK:', ', '.join(str(sp.nsimplify(e) if not e.free_symbols else sp.simplify(e)) for e in map(sp.sympify, expected)))
# ---
M, g, r, h = sp.symbols('M g r h', positive=True)
F = sp.Symbol('F', positive=True)
# 段の角 K を原点に：中心 O = (-s, r-h)、s = 水平距離
s = sp.sqrt(r**2 - (r - h)**2)
O = sp.Matrix([-s, r - h]); Tp = sp.Matrix([-s, 2 * r - h])
def tq(p, Fv): return p[0] * Fv[1] - p[1] * Fv[0]
grav = tq(O, sp.Matrix([0, -M * g]))
# 中心を K から遠ざかる向き？ 段は右にあるので右向き（+x）に押す
F1 = sp.solve(sp.Eq(tq(O, sp.Matrix([F, 0])) + grav, 0), F)[0]
check('MEC-RIG-205', '1', F1)
dd = sp.Symbol('dd', positive=True)   # dd = r - h > 0
Rk = sp.simplify(sp.sqrt(F1**2 + (M * g)**2).subs(h, r - dd)).subs(dd, r - h)
check('MEC-RIG-205', '2', Rk)
F2 = sp.solve(sp.Eq(tq(Tp, sp.Matrix([F, 0])) + grav, 0), F)[0]
check('MEC-RIG-205', '3', F2)
# 向きを変えて最小化：力の向き φ
phi = sp.Symbol('phi')
Fphi = sp.solve(sp.Eq(tq(Tp, sp.Matrix([F * sp.cos(phi), F * sp.sin(phi)])) + grav, 0), F)[0]
import numpy as np
for hv in (0.2, 0.5, 0.8):
    sub = {M: 1, g: 1, r: 1, h: hv}
    fn = sp.lambdify(phi, Fphi.subs(sub))
    ph = np.linspace(-1.5, 1.5, 300001); vals = fn(ph); vals = vals[vals > 0]
    assert abs(vals.min() - float(sp.sqrt(hv / 2))) < 1e-6
check('MEC-RIG-205', '4', M * g * sp.sqrt(h / (2 * r)))
