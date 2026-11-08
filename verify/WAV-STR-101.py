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
L, M, g_, rho = sp.symbols('L M g rho', positive=True)
v = sp.sqrt(M * g_ / rho)
f = 3 * v / (2 * L)                    # 3 loops: L = 3 λ/2
same('1', f, L=L, M=M, g=g_, rho=rho)
M2 = sp.symbols('M2', positive=True)
sol = sp.solve(sp.Eq(2 * sp.sqrt(M2 * g_ / rho) / (2 * L), f), M2)[0]
same('2', sol, M=M)
nums = {L: sp.Rational(12, 10), rho: sp.Rational(49, 10**5), M: sp.Rational(1, 2), g_: sp.Rational(98, 10)}
vv = v.subs(nums); assert vv == 100
same('3', vv)
ff = f.subs(nums); assert ff == 125
same('4', ff)
lam_s = sp.Rational(340) / ff
same('5', lam_s / 4)
# (6): frequency fixed by tuning fork; with M2, still 2 loops at the same f
f6 = 2 * sp.sqrt(sol * g_ / rho) / (2 * L)
assert sp.simplify(f6 - f) == 0
same('6', sp.Rational(340) / f6.subs(nums) / 4)
