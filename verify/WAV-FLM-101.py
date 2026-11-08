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
n, d, th, lam = sp.symbols('n d theta lam', positive=True)
# explicit geometric optical path difference: ray 2 inside film (2 d / cos r)·n minus ray 1 in air 2 d tan r sin θ
r = sp.asin(sp.sin(th) / n)
opd = 2 * d / sp.cos(r) * n - 2 * d * sp.tan(r) * sp.sin(th)
same('1', opd, d=d, n=n, theta=th)
# (2): both reflections low→high: phase π each → strong when OPD = mλ (checked in tex)
assert 'm\\lambda' in P['answer']['parts'][1]['tex']
# (3) normal incidence, minimum reflection: 2 n d = λ/2
nv, lv = 1.38, 552.0
d3 = lv / (4 * nv); same('3', round(d3, 9))
# phasor check: reflectivity at that thickness with Fresnel-like signs (both negative)
def refl(nf, ns, dd, lamb, theta=0.0):
    rr = math.asin(math.sin(theta) / nf)
    delta = 2 * math.pi / lamb * 2 * nf * dd * math.cos(rr)
    r1 = (1 - nf) / (1 + nf); r2 = (nf - ns) / (nf + ns)
    return abs(r1 + r2 * np.exp(-1j * delta))
ds = np.linspace(1, 300, 299001)
R = [refl(nv, 1.5, x, lv) for x in ds]
assert abs(ds[int(np.argmin(R))] - d3) < 0.01
# (4) θ = 60°, same order: λ' = 4 d sqrt(n²-sin²θ)
l4 = 4 * d3 * math.sqrt(nv**2 - math.sin(math.radians(60))**2)
same('4', num_sf(l4, 2))
ls = np.linspace(300, 700, 400001)
R4 = [refl(nv, 1.5, d3, x, math.radians(60)) for x in ls]
assert abs(ls[int(np.argmin(R4))] - l4) < 0.01
# (5) n = 2.0 > 1.5, maximum reflection
d5 = lv / (4 * 2.0)
same('5', num_sf(d5, 2))
R5 = [refl(2.0, 1.5, x, lv) for x in ds]
assert abs(ds[int(np.argmax(R5))] - d5) < 0.01
