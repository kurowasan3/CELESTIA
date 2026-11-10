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
v, V, m, M, dV = sp.symbols('v V m M dV', positive=True)
# 木星に対する速度は向きが逆転，大きさ不変
v1 = -(-v - V) + V
check('MEC-GRV-508', '2', v1)
g3 = grading('MEC-GRV-508', '3')
v1n = float(v1.subs({v: 10, V: 13}))
assert v1n == 36 and float(g3[0]) == 36
assert float(g3[1]) == sig((v1n / 10) ** 2), (v1n / 10) ** 2
dVs = sp.solve(sp.Eq(m * (-v) + M * V, m * v1 + M * (V - dV)), dV)[0]
check('MEC-GRV-508', '4', dVs)
# 第1次で一致
lossJ = sp.expand(M * V**2 / 2 - M * (V - dVs)**2 / 2)
gain = sp.expand(m * v1**2 / 2 - m * v**2 / 2)
assert sp.simplify(sp.limit((lossJ - gain) / m, M, sp.oo)) == 0
u = -(v - V) + V
check('MEC-GRV-508', '5', u)
import random
for _ in range(50):
    VV = random.uniform(1, 5); vv = VV + random.uniform(0.01, 10)
    assert abs(2 * VV - vv) < vv
print('MEC-GRV-508 OK')
