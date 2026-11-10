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
V = 340
v = sp.symbols('v', positive=True); f0 = sp.symbols('f0', positive=True)
s = sp.solve([sp.Eq(f0 * V / (V - v), 720), sp.Eq(f0 * V / (V + v), 640)], [v, f0], dict=True)[0]
same('1', s[v]); same('2', sf(float(s[f0]), 3))
# f は単調減少（近づく間も下がる）: 数値で確認
d = 34.0; vv = float(s[v]); F0 = float(s[f0])
import numpy as np
tau = np.linspace(-20, 20, 40001); r = np.hypot(vv * tau, d); t = tau + r / V
fobs = F0 * V / (V + vv * vv * tau / r)
assert np.all(np.diff(t) > 0) and np.all(np.diff(fobs) < 0); print('  (3) ok: 単調減少 → （イ）'); assert TEX['3'].startswith('（イ）')
same('4', round(vv * d / V, 6))
