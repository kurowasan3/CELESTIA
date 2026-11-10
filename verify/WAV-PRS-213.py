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
L = 3.0
lam = 4 * L / 3; same('1', sp.Rational(4) * sp.Symbol('L') / 3, L=sp.Symbol('L'))
x = np.linspace(0, L, 3000001); k = 2 * np.pi / lam
u = np.sin(k * x)                       # 変位の振幅分布（閉じた端 x=0 が節）
assert abs(abs(u[-1]) - 1) < 1e-9       # 開いた端が腹
nodes = x[np.abs(u) < 1e-6]; nodes = sorted(set(np.round(nodes, 3)))
same_set('2', nodes, L=L)
p = np.cos(k * x)                        # 密度変化 ∝ -du/dx の振幅分布
pmax = sorted(set(np.round(x[np.abs(np.abs(p) - 1) < 1e-12], 3))); same_set('3', pmax, L=L)
pz = sorted(set(np.round(x[np.abs(p) < 1e-6], 3))); same_set('4', pz, L=L)
# (5) ある瞬間 u(L/3) が正の最大：u = sin(kx)；密度変化 = -du/dx
dens = lambda xx: -k * np.cos(k * xx)
assert dens(0) < 0 and dens(2 * L / 3) > 0; print('  (5) ok: x=0 疎, 2L/3 密')
same('6', 425 * 4 * 0.60 / 3, rel=1e-9)
