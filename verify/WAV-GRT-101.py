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
import numpy as np
d = 1.0e-3 / 500
same('1', round(d * 1e6, 9))
lam = 600e-9
s1 = lam / d; same('2', round(s1, 12))
mmax = max(m for m in range(0, 50) if m * lam / d < 1)
same('3', 2 * mmax + 1)
x1 = 1.0 * math.tan(math.asin(s1))
same('4', num_sf(x1, 2))
assert num_sf(1.0 * s1, 2) != num_sf(x1, 2)    # small-angle approximation would give a different 2-sig-fig answer
# (5) overlap of 2nd (400-700) and 3rd orders
assert 2 * 700 > 3 * 400
print('  (5) ok: overlap 1200–1400 nm')
