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
pre = P['figures'][0]['sim']['preset']; A, lam, f = pre['A'], pre['lambda'], pre['f']
t1 = P['figures'][1]['sim']['t']; assert P['figures'][0]['sim']['t'] == 0 and t1 == 0.5
x = np.linspace(0, 8, 80001)
y0 = A * np.sin(2*np.pi*(-x/lam)); y1 = A * np.sin(2*np.pi*(f*t1 - x/lam))
same('1', float(y0.max()), i=0, rel=1e-6); same('1', lam, i=1)
# shift s in (0, λ) for +x direction: y1(x) = y0(x - s)
cands = [s for s in np.arange(0.001, lam, 0.001) if np.max(np.abs(A*np.sin(2*np.pi*(-(x-s)/lam)) - y1)) < 1e-3]
s = float(np.mean(cands)); v = s / t1; same('2', round(v, 3))
T = lam / v; same('3', round(T, 4))
ts = np.linspace(0, 5, 500001); yP = A*np.sin(2*np.pi*(ts/T - 2.0/lam))
i = np.argmax(yP > A - 1e-9); same('4', round(float(ts[i]), 3))
nn = sp.symbols('n', integer=True, nonnegative=True)
same('5', (sp.Rational(1) + 4*nn) / sp.Rational(1, 2), n=nn)
# -x direction with shift < λ
c2 = [s for s in np.arange(0.001, lam, 0.001) if np.max(np.abs(A*np.sin(2*np.pi*(-(x+s)/lam)) - y1)) < 1e-3]
same_set('6', [round(v, 3), round(float(np.mean(c2)) / t1, 3)])
