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
G, M, r = sp.symbols('G M r', positive=True)
vP, vQ = sp.symbols('v_P v_Q', positive=True)
check('MEC-GRV-204', '1', sp.sqrt(G * M / r))
sol = sp.solve([sp.Eq(r * vP, 3 * r * vQ), sp.Eq(vP**2 / 2 - G * M / r, vQ**2 / 2 - G * M / (3 * r))], [vP, vQ], dict=True)[0]
check('MEC-GRV-204', '2', sol[vP])
check('MEC-GRV-204', '3', sol[vQ])
a = (r + 3 * r) / 2
check('MEC-GRV-204', '4', sp.pi * sp.sqrt(a**3 / (G * M)))
check('MEC-GRV-204', '5', sp.sqrt(G * M / (3 * r)) / sol[vQ])
# 数値積分：P から vP で出発して，半周期後に距離 3r・速さ vQ になる
import numpy as np
GMn, rn = 1.0, 1.0
vp = float(sol[vP].subs({G: 1, M: 1, r: 1}))
y = np.array([rn, 0, 0, vp]); dt = 1e-4; T = float(sp.pi * sp.sqrt(8)); n = int(round(T / dt)); dt = T / n
def f(s):
    d = np.hypot(s[0], s[1]); return np.array([s[2], s[3], -s[0] / d**3, -s[1] / d**3])
for _ in range(n):
    k1 = f(y); k2 = f(y + dt / 2 * k1); k3 = f(y + dt / 2 * k2); k4 = f(y + dt * k3); y = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
assert abs(np.hypot(y[0], y[1]) - 3) < 1e-6 and abs(np.hypot(y[2], y[3]) - float(sol[vQ].subs({G: 1, M: 1, r: 1}))) < 1e-6
# シミュレーションの preset vratio = 1.22 ≈ sqrt(3/2)
assert abs(1.22 - float(sp.sqrt(sp.Rational(3, 2)))) < 0.005
