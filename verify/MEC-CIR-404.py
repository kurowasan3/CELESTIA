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
import math


def sig(x, n):
    """有効数字 n 桁に丸める"""
    x = float(x)
    if x == 0:
        return 0.0
    return round(x, n - 1 - int(math.floor(math.log10(abs(x)))))


def check_round(pid, part, *pairs):
    """数値の答え：独立に求めた値 e を有効数字 n 桁に丸めたものが grading の値と一致する（pairs = (e, n), ...）"""
    got = grading(pid, part)
    assert len(got) == len(pairs), (pid, part, got)
    for gv, (e, n) in zip(got, pairs):
        want = sig(sp.N(e), n)
        assert abs(float(gv) - want) <= 1e-12 * max(1, abs(want)), f'{pid} ({part}): grading {gv} ≠ {want}（{sp.N(e)} を {n} 桁）'
    print(f'{pid} ({part}) OK:', ', '.join(f'{float(sp.N(e)):.6g}→{sig(sp.N(e), n)}' for e, n in pairs))
# ---
m, g, l, d = sp.symbols('m g l d', positive=True)
v0 = sp.sqrt(2 * g * l)
check('MEC-CIR-404', '1', v0)
r = l - d
check('MEC-CIR-404', '2', m * g + m * v0**2 / l, m * g + m * v0**2 / r)
# 最高点（最下点から 2r）で張力 >= 0： v_t^2 >= g r
D = sp.Symbol('D', positive=True)
d0 = sp.solve(sp.Eq(v0**2 - 4 * g * (l - D), g * (l - D)), D)[0]
check('MEC-CIR-404', '3', d0)
# d = l/2：角 α での張力
a = sp.Symbol('alpha')
rr = l / 2
v2 = v0**2 - 2 * g * rr * (1 + sp.sin(a))
T = m * v2 / rr - m * g * sp.sin(a)
sa = sp.solve(sp.Eq(sp.simplify(T.subs(sp.sin(a), sp.Symbol('S'))), 0), sp.Symbol('S'))[0]
vv = sp.sqrt(v2.subs(sp.sin(a), sa))
check('MEC-CIR-404', '4', sa, vv)
# たるんだ後の放物運動の最高点（最下点から）
ca = sp.sqrt(1 - sa**2)
y0 = rr * (1 + sa)
Hmax = y0 + (vv * ca)**2 / (2 * g)
check('MEC-CIR-404', '5', sp.simplify(l - Hmax))   # 2l/27
# 数値で：P のまわりの円運動を積分し，張力 0 の角を確かめる
import numpy as np
gg, ll = 9.8, 1.0
R_ = ll / 2; phi, om, dt = 0.0, np.sqrt(2 * gg * ll) / R_, 1e-6   # phi：最下点から測った角
while True:
    Tn = R_ * om**2 + gg * np.cos(phi)          # 単位質量あたりの張力（中心向き）
    if Tn <= 0: break
    om += -gg / R_ * np.sin(phi) * dt; phi += om * dt
alpha = phi - np.pi / 2                           # 水平線から上に測った角
assert abs(np.sin(alpha) - 2 / 3) < 1e-4
