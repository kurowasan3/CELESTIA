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
M, m, g, l, x, mu, theta = sp.symbols('M m g l x mu theta', positive=True)
N, R, f = sp.symbols('N R f')
# 力のつり合い＋下端 A のまわりのモーメント（ベクトルの外積で）
def tq(px, py, Fx, Fy): return px * Fy - py * Fx
c, s = sp.cos(theta), sp.sin(theta)
eqs = [sp.Eq(N - (M + m) * g, 0), sp.Eq(f - R, 0),
       sp.Eq(tq(l / 2 * c, l / 2 * s, 0, -M * g) + tq(x * c, x * s, 0, -m * g) + tq(l * c, l * s, -R, 0), 0)]
sol = sp.solve(eqs, [N, R, f], dict=True)[0]
check('MEC-RIG-204', '1', sol[N])
check('MEC-RIG-204', '2', sol[R], sol[f])
xmax = sp.solve(sp.Eq(sol[f], mu * sol[N]), x)[0]
check('MEC-RIG-204', '3', xmax)
# (4) x=l で f <= mu N  ⇔  tanθ >= (M+2m)/(2μ(M+m))
tcrit = sp.solve(sp.Eq(sol[f].subs(x, l), mu * sol[N]), theta)[0]
assert same(sp.tan(tcrit), (M + 2 * m) / (2 * mu * (M + m)))
