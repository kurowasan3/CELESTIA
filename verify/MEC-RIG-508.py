# --- 共通：problem.json の入力用の答え（grading）を sympy の式に直し，独立に導いた式と数値で照合する
import json, re, random, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor
ROOT = pathlib.Path(__file__).resolve().parent.parent
TR = standard_transformations + (implicit_multiplication_application, convert_xor)
GREEK = {'π': 'pi', 'θ': 'theta', 'μ': 'mu', 'ω': 'omega', 'φ': 'phi', 'α': 'alpha', 'ρ': 'rho'}


def raw(pid, part):
    f = next(ROOT.glob(f'content/problems/*/{pid}/problem.json'))
    p = json.loads(f.read_text(encoding='utf-8'))
    return next(x for x in p['answer']['parts'] if x['id'] == part)['grading']


def grading(pid, part):
    out = []
    for fld in re.split(r'\s;\s', raw(pid, part)):
        if '：' in fld:
            fld = fld.split('：', 1)[1]
        fld = re.sub(r'\[[^\]]*\]\s*$', '', fld).strip()
        fld = re.sub(r'(倍|°|枚)$', '', fld).strip()
        s = fld
        for k, v in GREEK.items():
            s = s.replace(k, v)
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
            if n in ('theta', 'phi', 'alpha'):
                vals[n] = random.uniform(0.2, 1.3)
            elif n in ('mu', 'mu_0'):
                vals[n] = random.uniform(0.05, 0.3)
            else:
                vals[n] = random.uniform(1.0, 3.0)
        sub = {s: vals[s.name] for s in allsyms}
        x, y = complex(a.subs(sub).evalf()), complex(b.subs(sub).evalf())
        if abs(x.imag) > 1e-12 or abs(y.imag) > 1e-12:
            continue
        ok += 1
        if abs(x - y) > 1e-9 * max(1, abs(x), abs(y)):
            return False
    return ok >= trials


def check(pid, part, *expected):
    got = grading(pid, part)
    assert len(got) == len(expected), (pid, part, got)
    for g, e in zip(got, expected):
        assert same(g, e), f'{pid} ({part}): grading {g} ≠ 導出 {sp.simplify(e)}'
    print(f'{pid} ({part}) OK:', ', '.join(str(sp.simplify(sp.sympify(e))) for e in expected))


def sig(x, n=2):
    """有効数字 n 桁に丸めた値"""
    x = float(x)
    if x == 0:
        return 0.0
    return float(f'{x:.{n - 1}e}')
# ---
R, r, th, T, m, g, mu, f, N = sp.symbols('R r theta T m g mu f N', positive=True)
# 接点 P を原点，中心 C=(0,R)。糸は軸の下側から，向き (cos th, sin th)。接点 = C + r(sin th, -cos th)
C = sp.Matrix([0, R]); Tg = C + r * sp.Matrix([sp.sin(th), -sp.cos(th)]); d = sp.Matrix([sp.cos(th), sp.sin(th)])
torqueP = sp.simplify(Tg[0] * T * d[1] - Tg[1] * T * d[0])
assert sp.simplify(torqueP - T * (r - R * sp.cos(th))) == 0
assert torqueP.subs(th, 0).subs({r: 1, R: 2, T: 1}) < 0            # θ=0：時計回り → 引く向き（+x）に転がる
c0 = sp.solve(sp.Eq(r - R * sp.Symbol('c'), 0), sp.Symbol('c'))[0]
check('MEC-RIG-508', '2', c0)
th0 = sp.acos(r / R)
fs = T * sp.cos(th0); Ns = m * g - T * sp.sin(th0)
assert sp.simplify(T * r - fs * R) == 0                            # 中心のまわりも 0
check('MEC-RIG-508', '3', fs, Ns)
T1 = sp.solve(sp.Eq(fs, mu * Ns), T)[0]
check('MEC-RIG-508', '4', T1)
T2 = sp.solve(sp.Eq(Ns, 0), T)[0]
assert same(sp.simplify(T2 - T1) > 0, True) if False else True
import random
for _ in range(100):
    v = {R: random.uniform(1, 3), mu: random.uniform(0.01, 3), m: 1, g: 9.8}
    v[r] = v[R] * random.uniform(0.01, 0.99)
    assert float(T1.subs(v)) < float(T2.subs(v))
assert torqueP.subs({th: 1.4, r: 1, R: 2, T: 1}) > 0                # θ>θ0：反時計回り → 遠ざかる
print('MEC-RIG-508 OK')
