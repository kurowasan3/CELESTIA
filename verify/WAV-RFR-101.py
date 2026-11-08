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
c, lam0, n = sp.Rational(3, 1) * 10**8, sp.Rational(600) * sp.Rational(1, 10**9), sp.Rational(3, 2)
f = c / lam0
same('1', f)
v = c / n; same('2', v)
lam_g = v / f                          # frequency unchanged across boundary
same('3', lam_g * 10**9)
same('4', f)
assert sp.Rational(1, 3) * 10**15 == v / lam0  # the listed wrong value (keeping λ, slowing v)
s2 = sp.sin(sp.pi / 6) / n
same('5', s2)
pre = P['simulation']['preset']
assert abs(math.degrees(math.asin(math.sin(math.radians(pre['theta1'])) * pre['n1'] / pre['n2'])) - math.degrees(math.asin(1/3))) < 1e-9
