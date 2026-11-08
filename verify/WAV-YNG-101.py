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
L, lam, d, n, t, s, l = sp.symbols('L lam d n t s l', positive=True)
same('1', L * lam / d, L=L, lam=lam, d=d)
x = sp.symbols('x', real=True)
# exact geometry numerically: S1 at y=+d/2, S2 at y=-d/2; screen at distance L; x upward
def zero_order(dv, Lv, lv, extra1=0.0, s0=0.0):
    def opd(xx):
        p1 = math.hypot(lv, s0 - dv / 2) + math.hypot(Lv, xx - dv / 2) + extra1
        p2 = math.hypot(lv, s0 + dv / 2) + math.hypot(Lv, xx + dv / 2)
        return p1 - p2
    a, b = -0.1, 0.1
    for _ in range(200):
        m = (a + b) / 2
        if opd(a) * opd(m) <= 0: b = m
        else: a = m
    return (a + b) / 2
dv, Lv, lv = 0.3e-3, 1.5, 0.20
nv, tv, lamv = 1.5, 6.0e-6, 600e-9
x2 = zero_order(dv, Lv, lv, extra1=(nv - 1) * tv)
want2 = g('2', n=n, t=t, L=L, d=d).subs({n: nv, t: tv, L: Lv, d: dv})
assert abs(x2 - float(want2)) / abs(x2) < 1e-3, (x2, want2)
same('2', (n - 1) * t * L / d, n=n, t=t, L=L, d=d)
same('3', sp.nsimplify((nv - 1) * tv / lamv))
sv = 1.0e-3
x4 = zero_order(dv, Lv, lv, s0=sv)
want4 = g('4', s=s, L=L, l=l).subs({s: sv, L: Lv, l: lv})
assert abs(x4 - float(want4)) / abs(x4) < 1e-3, (x4, want4)
same('4', -s * L / l, s=s, L=L, l=l)
s5 = (nv - 1) * tv * lv / dv
assert abs(zero_order(dv, Lv, lv, extra1=(nv - 1) * tv, s0=s5)) < 1e-6
same('5', round(s5 * 1e3, 9))
