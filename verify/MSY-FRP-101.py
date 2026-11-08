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
v_0, mu, g, L, H, D = sp.symbols('v_0 mu g L H D', positive=True)
m = sp.Symbol('m', positive=True)
v1 = sp.sqrt(sp.solve(sp.Eq(m * sp.Symbol('w')**2 / 2 - m * v_0**2 / 2, -mu * m * g * L), sp.Symbol('w'))[1]**2)
check('MSY-FRP-101', '1', v1)
t = sp.sqrt(2 * H / g)
check('MSY-FRP-101', '2', t)
check('MSY-FRP-101', '3', v1 * t)
mus = sp.solve(sp.Eq(D, v1 * t), mu)[0]
check('MSY-FRP-101', '4', mus)
val = mus.subs({v_0: 3, L: sp.Rational(1, 2), H: sp.Rational(4, 5), g: sp.Rational(98, 10), D: sp.Rational(4, 5)})
print('  μ =', float(val))
assert round(float(val), 2) == 0.52
check('MSY-FRP-101', '5', sp.Rational(52, 100))
