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
# 動摩擦のあるばね振り子：区分ごとに運動方程式を数値積分し，式とも照合
m, k, g, A, mu, mu_0 = sp.symbols('m k g A mu mu_0', positive=True)
x1 = sp.solve(sp.Eq(-k * sp.Symbol('x') + mu * m * g, 0), sp.Symbol('x'))[0]
check('MEC-SHM-203', '1', x1)
check('MEC-SHM-203', '2', 2 * x1 - A)
check('MEC-SHM-203', '3', sp.pi * sp.sqrt(m / k))
# (4) A <= (mu_0 + 2 mu) m g / k  ⇔  k|x2| <= mu_0 m g
Acrit = sp.solve(sp.Eq(k * (A - 2 * x1), mu_0 * m * g), A)[0]
assert sp.simplify(Acrit - (mu_0 + 2 * mu) * m * g / k) == 0
# (5) 数値シミュレーション
import numpy as np
mm, kk, mk, ms, AA, gg = 0.5, 49.0, 0.20, 0.30, 0.10, 9.8
X, V, t, dt = AA, 0.0, 0.0, 2e-6
moving = True
while True:
    if V == 0.0:
        if abs(kk * X) <= ms * mm * gg: break
        sgn = -np.sign(X)          # 動き出す向き
    else:
        sgn = np.sign(V)
    a_ = (-kk * X - sgn * mk * mm * gg) / mm
    Vn = V + a_ * dt
    if V != 0.0 and np.sign(Vn) != np.sign(V):
        V = 0.0; X += 0; t += dt; continue
    V = Vn; X += V * dt; t += dt
assert abs(X - 0.02) < 1e-4, X
assert abs(t - float(sp.sqrt(2) * sp.pi / 7)) < 1e-3, t
check('MEC-SHM-203', '5', sp.Rational(2, 100), 2 * sp.pi * sp.sqrt(sp.Rational(1, 2) / 49))
