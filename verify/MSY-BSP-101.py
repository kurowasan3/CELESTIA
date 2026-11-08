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
m, M, k, v = sp.symbols('m M k v', positive=True)
u = m * v / (m + M)
check('MSY-BSP-101', '1', u)
VG = m * v / (m + 2 * M)
check('MSY-BSP-101', '2', VG)
x = sp.Symbol('x', positive=True)
xm = sp.solve(sp.Eq((m + M) * u**2 / 2, (m + 2 * M) * VG**2 / 2 + k * x**2 / 2), x)[0]
check('MSY-BSP-101', '3', xm)
mu_r = (m + M) * M / (m + 2 * M)
check('MSY-BSP-101', '4', sp.pi / 2 * sp.sqrt(mu_r / k))
# 2体の運動を解析的に：相対座標 s = xB - xA（ばねの伸び），s'' = -k s/mu_r
t = sp.Symbol('t', nonnegative=True)
w = sp.sqrt(k / mu_r)
s = -(u / w) * sp.sin(w * t)              # s(0)=0, s'(0) = 0 - u
vA = VG - M / (m + 2 * M) * sp.diff(s, t)  # xA = XG - M/(m+2M) s
vB = VG + (m + M) / (m + 2 * M) * sp.diff(s, t)
assert same(vA.subs(t, 0), u) and same(vB.subs(t, 0), 0)
check('MSY-BSP-101', '5', vB.subs(t, sp.pi / w), vA.subs(t, sp.pi / w))
# 数値：極値であることを確認
import numpy as np
sub = {m: 0.05, M: 1.0, k: 50.0, v: 100.0}
fa = sp.lambdify(t, vA.subs(sub)); fb = sp.lambdify(t, vB.subs(sub))
ts = np.linspace(0, 5, 200001)
assert abs(fa(ts).min() - float(vA.subs(t, sp.pi / w).subs(sub))) < 1e-6 and fa(ts).min() > 0
assert abs(fb(ts).max() - float(vB.subs(t, sp.pi / w).subs(sub))) < 1e-6
