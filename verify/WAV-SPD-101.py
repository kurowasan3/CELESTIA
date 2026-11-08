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
V = 340
same('1', sp.Rational(V, 425))
same('2', V * 2.0 / 2)
same('3', 425)                                   # 振動数は音源で決まる
same('4', num_sf(346 / 425, 2))
