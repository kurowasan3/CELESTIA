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
# (1) on segment, in phase: |l1-l2| = |2x| = mλ
d, lam, x, delta = sp.symbols('d lam x delta', positive=True)
# (2),(3): count hyperbolae |Δ| = c λ with |Δ| < d
def count(dd, offset):
    return sum(1 for m in range(-20, 21) if abs(m + offset) < dd)
same('2', count(2.2, 0)); same('3', count(2.2, 0.5))
# brute-force check of (2): number of amplitude maxima along a segment parallel... use the segment itself
xs = np.linspace(-1.1 + 1e-6, 1.1 - 1e-6, 200001)            # λ=1, d=2.2
for ph, want in ((0, 5), (np.pi, 4)):
    l1, l2 = 1.1 + xs, 1.1 - xs
    amp = np.abs(np.exp(-1j * 2 * np.pi * l1) + np.exp(-1j * (2 * np.pi * l2 + ph)))
    peaks = np.sum((amp[1:-1] > amp[:-2]) & (amp[1:-1] >= amp[2:]) & (amp[1:-1] > 1.99))
    assert peaks == want, (ph, peaks)
print('  (2)(3) brute force ok')
# (4) S2 lags by δ: phase S1 = ωt - k l1, S2 = ωt - δ - k l2
k = 2 * sp.pi / lam
l1, l2 = d / 2 + x, d / 2 - x
cond = sp.solve(sp.Eq((-k * l1) - (-delta - k * l2), 0), x)[0]
same('4', cond, delta=delta, lam=lam)
assert cond.subs({delta: 1, lam: 1}) > 0          # toward S2 (x>0)
# numeric amplitude-maximum position check
lv, dl = 2.0, math.pi / 2
xs = np.linspace(-0.6, 0.6, 1200001)
amp = np.abs(np.exp(-1j * 2 * np.pi / lv * (2.2 + xs)) + np.exp(-1j * (2 * np.pi / lv * (2.2 - xs) + dl)))
xm = xs[np.argmax(amp)]
same('5', round(float(xm), 4))
assert abs(xm - 0.25) < 1e-5
pre = P['simulation']['preset']; assert abs(pre['d'] / pre['lambda'] - 2.2) < 1e-12 and pre['phase'] == 0
