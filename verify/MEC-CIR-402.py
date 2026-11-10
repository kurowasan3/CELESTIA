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
g = sp.Rational(49, 5); m = sp.Rational(1, 5); r = sp.Rational(1, 2); h = sp.Rational(5, 2)
vB = sp.sqrt(2 * g * h)
check('MEC-CIR-402', '1', vB)                     # 7
NB = m * g + m * vB**2 / r
check_round('MEC-CIR-402', '2', (NB, 2))          # 21.56 -> 22
vC = sp.sqrt(2 * g * (h - 2 * r))
NC = m * vC**2 / r - m * g
check_round('MEC-CIR-402', '3', (vC, 2), (NC, 2))  # 5.42 -> 5.4, 9.8
assert NC == sp.Rational(49, 5)
H = sp.Symbol('H', positive=True)
h0 = sp.solve(sp.Eq(2 * g * (H - 2 * r), g * r), H)[0]
check('MEC-CIR-402', '4', h0)                     # 1.25
