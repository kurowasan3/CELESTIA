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
# 鉛直ばね上の板と小物体：x はつり合いからの上向きの変位
M, m, k, g, a = sp.symbols('M m k g a', positive=True)
x = sp.Symbol('x')
d0 = sp.solve(sp.Eq(k * sp.Symbol('d'), (M + m) * g), sp.Symbol('d'))[0]
check('MEC-SHM-202', '1', d0)
acc = (k * (d0 - x) - (M + m) * g) / (M + m)        # = -k x/(M+m)
w2 = sp.simplify(-acc / x)
check('MEC-SHM-202', '2', 2 * sp.pi / sp.sqrt(w2))
N = m * (g + acc)
check('MEC-SHM-202', '3', N.subs(x, -a))
# 離れない条件：最高点 x=a で N>=0  →  a <= (M+m)g/k
amax = sp.solve(sp.Eq(N.subs(x, a), 0), a)[0]
assert sp.simplify(amax - (M + m) * g / k) == 0
# N=0 となる位置はばねが自然長（縮み d0 - x = 0）
xs = sp.solve(sp.Eq(N, 0), x)[0]
assert sp.simplify(d0 - xs) == 0
v = sp.sqrt(w2 * (a**2 - xs**2))
check('MEC-SHM-202', '5', v)
# 数値積分でも確認（板と小物体を一体で）
import numpy as np
MM, mm, kk, gg, aa = 1.0, 0.5, 100.0, 9.8, 0.25
d = (MM + mm) * gg / kk; X, Vv, dt = -aa, 0.0, 1e-6
while True:
    acc_ = -kk * X / (MM + mm)
    if gg + acc_ <= 0: break
    Vv += acc_ * dt; X += Vv * dt
assert abs(X - d) < 1e-3 and abs(Vv - float(v.subs({M: MM, m: mm, k: kk, g: gg, a: aa}))) < 1e-3
