import json, re, math, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor
ROOT = Path(__file__).resolve().parent.parent
PID = Path(__file__).stem
P = json.loads((ROOT / 'content/problems/WAV' / PID / 'problem.json').read_text(encoding='utf-8'))
G = {a['id']: a.get('grading') for a in P['answer']['parts']}
GREEK = {'π': ' pi ', 'λ': ' lam ', 'ρ': ' rho ', 'δ': ' delta ', 'θ': ' theta '}
def g(pid, **syms):
    """grading 文字列を sympy 式に（ラベル・単位・倍/本/個/回 を外す）"""
    s = G[pid]
    assert s, f'({pid}) に grading がない'
    s = re.sub(r'^[^：]*：', '', s)
    s = re.sub(r'\s*\[[^\]]*\]\s*$', '', s)
    s = re.sub(r'(倍|本|個|回)$', '', s)
    for k, v in GREEK.items():
        s = s.replace(k, v)
    loc = {'pi': sp.pi, 'sqrt': sp.sqrt, 'sin': sp.sin, 'cos': sp.cos}
    loc.update(syms)
    return parse_expr(s, local_dict=loc, transformations=standard_transformations + (implicit_multiplication, convert_xor))
def same(pid, want, rel=1e-9, **syms):
    got = g(pid, **syms)
    free = sorted(sp.sympify(got).free_symbols | sp.sympify(want).free_symbols, key=str)
    for _ in range(6):
        sub = {s: sp.Rational(random.randint(11, 97), random.randint(5, 13)) for s in free}
        a, b = complex(sp.N(sp.sympify(got).subs(sub))), complex(sp.N(sp.sympify(want).subs(sub)))
        assert abs(a - b) <= rel * max(1, abs(a), abs(b)), f'({pid}) grading {G[pid]!r} = {a} != {b}'
    print(f'  ({pid}) ok: {G[pid]}')
def num_sf(x, n):
    """有効数字 n 桁に丸める"""
    if x == 0: return 0.0
    return float(f'{x:.{n - 1}e}')
print(PID)
x, t = sp.symbols('x t', real=True)
y = sp.Rational(1, 5) * sp.sin(sp.pi * (4 * t - sp.Rational(1, 2) * x))
A = sp.Rational(1, 5)
omega = sp.diff(4 * sp.pi * t, t); k = sp.pi / 2
f = omega / (2 * sp.pi); lam = 2 * sp.pi / k
same('1', A); same('2', f); same('3', lam)
v = omega / k
same('4', v)
# direction: phase const → x = 8 t + c  (dx/dt > 0)
assert sp.solve(sp.Eq(4 * t - x / 2, 0), x)[0].diff(t) > 0
y3 = y.subs(x, 3)
assert sp.simplify(y3 - sp.Rational(1, 5) * sp.cos(4 * sp.pi * t)) == 0
same('5', sp.Rational(1, 5) * sp.cos(4 * sp.pi * t), t=t)
y0 = y.subs(t, 0)
assert sp.simplify(y0 + sp.Rational(1, 5) * sp.sin(sp.pi * x / 2)) == 0
# answer figure preset matches the equation: A sin 2π(f t − x/λ)
pre = P['figures'][0]['sim']['preset']
assert (pre['A'], pre['f'], pre['lambda'], pre['dir']) == (0.2, 2, 4, 1)
print('  (6) ok: y = -0.2 sin(πx/2), preset matches')
