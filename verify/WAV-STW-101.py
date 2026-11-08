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
A, T, lam, L, x, t = sp.symbols('A T lam L x t', positive=True)
y1 = A * sp.sin(2 * sp.pi * (t / T - x / lam))
y2 = g('1', A=A, T=T, lam=lam, L=L, x=x, t=t)
# reflected wave travels toward -x: function of t/T + x/λ only
assert sp.simplify(sp.diff(y2, t) * T - sp.diff(y2, x) * lam) == 0
# fixed end
assert sp.simplify(sp.expand_trig((y1 + y2).subs(x, L))) == 0
print('  (1) ok: y2 left-moving and y1+y2=0 at x=L')
ys = g('2', A=A, T=T, lam=lam, L=L, x=x, t=t)
for _ in range(20):
    sub = {A: random.random(), T: random.random() + .1, lam: random.random() + .1, L: random.random() * 3, x: random.random() * 3, t: random.random() * 5}
    assert abs(float((y1 + y2 - ys).subs(sub))) < 1e-12
print('  (2) ok: product form equals y1+y2')
# (3) nodes x = L - mλ/2
for m in range(5):
    xm = L - m * lam / 2
    assert sp.simplify(sp.sin(2 * sp.pi * (L - xm) / lam)) == 0
# (4) max speed at antinode
xa = L - lam / 4
vmax = sp.Abs(sp.diff(ys, t).subs(x, xa))
same('4', 4 * sp.pi * A / T, A=A, T=T)
import numpy as np
vfun = sp.lambdify(t, sp.diff(ys, t).subs(x, xa).subs({A: 0.3, T: 1.7, lam: 0.9, L: 2.3}), 'numpy')
assert abs(np.abs(vfun(np.linspace(0, 1.7, 200001))).max() - 4 * math.pi * 0.3 / 1.7) < 1e-6
# (5) count nodes on 0<=x<=L numerically
Lv, lv = 3.0, 0.8
nodes = [Lv - m * lv / 2 for m in range(100) if Lv - m * lv / 2 >= -1e-12]
same('5', len(nodes))
# (6) amplitude at x=0
amp_fixed = abs(2 * math.sin(2 * math.pi * Lv / lv))
amp_free = abs(2 * math.cos(2 * math.pi * Lv / lv))     # free end: y2 sign flipped → 2A cos(...)
assert abs(amp_fixed - 2) < 1e-9 and amp_free < 1e-9
assert '固定端のとき腹，自由端のとき節' in P['answer']['parts'][5]['tex']
print('  (6) ok: fixed→antinode, free→node at x=0')
