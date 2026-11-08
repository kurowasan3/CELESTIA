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
# y-t graph of x=0 from the sim preset (wave.traveling: y = A sin 2π(f t − x/λ))
import numpy as np
pre = P['figures'][0]['sim']['preset']
assert pre == P['simulation']['preset'] and pre['xP'] == 0 and P['figures'][0]['plot'] == {'x': 't', 'y': ['yP']}
t = np.linspace(0, 6, 600001)
y = pre['A'] * np.sin(2 * np.pi * (pre['f'] * t - pre['xP'] / pre['lambda']))
A = y.max(); same('1', round(A, 6))
up = t[1:][(y[:-1] < 0) & (y[1:] >= 0)]          # upward zero crossings
T = float(np.diff(up).mean()); same('2', round(T, 6))
same('3', round(1 / T, 6))
v = 2.0
lam = v * T; same('4', round(lam, 6))
assert abs(lam - pre['lambda']) < 1e-9, 'preset λ と v T が食い違う'
i = np.searchsorted(t, 1.0); assert abs(y[i]) < 1e-6 and y[i + 10] < 0   # (5) 負の向き
assert '負' in P['answer']['parts'][4]['tex']
print('  (5) ok: 負の向き')
