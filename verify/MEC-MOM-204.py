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
# はね返り続けるボール：等比級数の和と数値シミュレーション
m, g, h, e = sp.symbols('m g h e', positive=True)
v0 = sp.sqrt(2 * g * h); t0 = sp.sqrt(2 * h / g)
check('MEC-MOM-204', '1', e * v0)
check('MEC-MOM-204', '2', 2 * e * v0 / g)
n = sp.Symbol('n', integer=True, nonnegative=True)
T = t0 + sp.summation(2 * e**(n + 1) * v0 / g, (n, 0, sp.oo)).subs(sp.Abs(e), e)
T = sp.simplify(sp.piecewise_fold(T).args[0][0]) if isinstance(T, sp.Piecewise) else T
check('MEC-MOM-204', '3', T)
S = h + sp.summation(2 * e**(2 * n + 2) * h, (n, 0, sp.oo))
S = S.args[0][0] if isinstance(S, sp.Piecewise) else S
S = sp.piecewise_fold(S); S = S.args[0][0] if isinstance(S, sp.Piecewise) else S
check('MEC-MOM-204', '4', S)
J = sp.summation(m * (1 + e) * e**n * v0, (n, 0, sp.oo))
J = sp.piecewise_fold(J); J = J.args[0][0] if isinstance(J, sp.Piecewise) else J
assert same(J, m * g * T)          # 衝突ごとの力積の和 = mgT
check('MEC-MOM-204', '5', m * g * T)
# 数値：衝突を 400 回たどる
gg, hh, ee = 9.8, 1.0, 0.6
v = (2 * gg * hh) ** 0.5; tot = (2 * hh / gg) ** 0.5; dist = hh
for _ in range(400):
    v *= ee; tot += 2 * v / gg; dist += v * v / gg
sub = {g: gg, h: hh, e: ee}
assert abs(tot - float(T.subs(sub))) < 1e-9 and abs(dist - float(S.subs(sub))) < 1e-9
