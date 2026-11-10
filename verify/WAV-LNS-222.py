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
data = [(15, 29.5), (20, 20.3), (25, 16.5), (30, 15.2), (40, 13.4)]
X = np.array([1 / a for a, b in data]); Y = np.array([1 / b for a, b in data])
# 傾き -1 の直線を当てはめた切片
c = float(np.mean(Y + X)); print('  intercept (slope -1 fit):', c)
slope, ic = np.polyfit(X, Y, 1); print('  free fit slope, intercept:', slope, ic)
assert abs(c - 0.10) < 0.002 and abs(slope + 1) < 0.05
same('2', round(1 / round(c, 2), 6))
b = sp.Rational(1) / (sp.Rational(1, 10) - sp.Rational(1, 15)); same('3', b / 15)
s, a, f = sp.symbols('s a f', positive=True)
disc = sp.discriminant(a**2 - s * a + s * f, a); assert sp.simplify(disc - (s**2 - 4 * s * f)) == 0
print('  (5) ok: s >= 4f =', 4 * 10)
spec = json.loads((ROOT / 'content/problems/WAV' / PID / 'q-graph.fig.json').read_text(encoding='utf-8'))
pts = [it['at'] for it in spec['items'] if it['t'] == 'circle']
assert len(pts) == 5; print('  figure points:', len(pts))
