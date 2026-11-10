import json, re, math, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor
ROOT = Path(__file__).resolve().parent.parent
PID = Path(__file__).stem
P = json.loads((ROOT / 'content/problems/WAV' / PID / 'problem.json').read_text(encoding='utf-8'))
G = {a['id']: a.get('grading') for a in P['answer']['parts']}
TEX = {a['id']: a['tex'] for a in P['answer']['parts']}
GREEK = {'π': 'pi', 'λ': 'lam', 'ρ': 'rho', 'δ': 'delta', 'θ': 'theta', 'Δ': 'Delta'}
def field(pid, i=0):
    s = G[pid]
    assert s, f'({pid}) に grading がない'
    s = [x.strip() for x in s.split(';')][i]
    s = re.sub(r'^[^：]*：', '', s)
    s = re.sub(r'\s*\[[^\]]*\]\s*$', '', s)
    s = re.sub(r'(倍|本|個|回|°)$', '', s)
    return s.strip()
def parse(s, **syms):
    for k, v in GREEK.items():
        s = re.sub(k + r'(?=_)', v, s)
        s = s.replace(k, f' {v} ')
    loc = {'pi': sp.pi, 'sqrt': sp.sqrt, 'sin': sp.sin, 'cos': sp.cos}
    loc.update(syms)
    return parse_expr(s, local_dict=loc, transformations=standard_transformations + (implicit_multiplication, convert_xor))
def close(a, b, rel):
    return abs(a - b) <= rel * max(1, abs(a), abs(b))
def same(pid, want, i=0, rel=1e-9, **syms):
    got = parse(field(pid, i), **syms)
    want = sp.sympify(want)
    free = sorted(sp.sympify(got).free_symbols | want.free_symbols, key=str)
    done = 0
    while done < 6:
        sub = {s: sp.Rational(random.randint(11, 97), random.randint(5, 13)) for s in free}
        ea, eb = sp.N(got.subs(sub)), sp.N(want.subs(sub))
        if not (ea.is_finite and eb.is_finite):
            continue   # 分母が 0 になる代入は除く
        a, b = complex(ea), complex(eb)
        assert close(a, b, rel), f'({pid}) grading {G[pid]!r} = {a} != {b}'
        done += 1
    print(f'  ({pid}) ok: {G[pid]}')
def same_set(pid, wants, i=0, rel=1e-9, **syms):
    s = field(pid, i)
    assert s.startswith('{') and s.endswith('}')
    got = sorted(complex(sp.N(sp.sympify(parse(x, **syms)))).real for x in s[1:-1].split(','))
    assert len(got) == len(wants) and all(close(a, b, rel) for a, b in zip(got, sorted(wants))), f'({pid}) {got} != {wants}'
    print(f'  ({pid}) ok: {G[pid]}')
def sf(x, n):
    """有効数字 n 桁"""
    return float(f'{x:.{n - 1}e}')
print(PID)
import numpy as np
pre = P['figures'][0]['sim']['preset']
assert P['figures'][0]['sim']['t'] == 0 and pre['dir'] == 1
A, lam, f = pre['A'], pre['lambda'], pre['f']
x = np.linspace(0, 8, 800001)
y = A * np.sin(2 * np.pi * (f * 0 - x / lam))
# wavelength: distance between successive downward... use crests
cr = x[1:-1][(y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:])]
same('1', round(float(np.diff(cr).mean()), 4))
# density excess ∝ -dy/dx
dens = -np.gradient(y, x)
inner = (x > 0.01) & (x < 7.99)
dmax = x[inner][np.argmax(dens[inner])]; same('2', round(dmax, 3))
xs = x[inner]; dd = dens[inner]
mins = xs[1:-1][(dd[1:-1] < dd[:-2]) & (dd[1:-1] <= dd[2:])]
same_set('3', sorted(set(round(m, 2) for m in mins)))
vel = A * 2 * np.pi * f * np.cos(2 * np.pi * (-x / lam))
zero = x[1:][(np.sign(vel[:-1]) != np.sign(vel[1:]))]
same_set('4', [round(z, 2) for z in zero if 0 < z < 4])
i = np.searchsorted(x, 2.0); assert vel[i] < 0 and '負' in TEX['5']; print('  (5) ok: 負の向き')
