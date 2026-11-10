# 検算スクリプト（problem.json の入力用の答えを sympy で独立に導いた値と照合する）
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ELM-LOR-402'
PROB = json.loads((Path(__file__).resolve().parent.parent / 'content' / 'problems' / PID[:3] / PID / 'problem.json').read_text(encoding='utf-8'))
GRADING = {a['id']: a.get('grading') for a in PROB['answer']['parts']}
GREEK = {'μ': 'mu', 'θ': 'theta', 'ω': 'omega', 'Φ': 'Phi', 'τ': 'tau', 'Δ': 'Delta', 'φ': 'phi'}
FUNCS = {'abs', 'sqrt', 'sin', 'cos', 'tan', 'log', 'exp', 'pi'}
TR = standard_transformations + (implicit_multiplication, convert_xor)


def S(names):
    """記号（名前は入力用の答えと同じ：μ_0 → mu_0）"""
    return sp.symbols(names, positive=True)


def parse(src):
    s = src.strip()
    for k, v in GREEK.items():
        s = s.replace(k, ' ' + v + ' ' if v != 'mu' else v)
    s = s.replace('π', ' pi ').replace('×', '*')
    s = re.sub(r'\s+_', '_', s)
    s = re.sub(r'(\d)\s*([A-Za-z(])', r'\1*\2', s)
    names = set(re.findall(r'[A-Za-z][A-Za-z0-9_]*', s)) - FUNCS
    loc = {n: sp.Symbol(n, positive=True) for n in names}
    loc.update({'abs': sp.Abs, 'sqrt': sp.sqrt, 'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'log': sp.log, 'exp': sp.exp, 'pi': sp.pi})
    return parse_expr(s, local_dict=loc, transformations=TR)


def field(pid, k=0):
    g = GRADING[pid]
    assert g, f'({pid}) に入力用の答えがありません'
    f = [x.strip() for x in g.split(' ; ')][k]
    if '：' in f:
        f = f.split('：', 1)[1]
    f = re.sub(r'\s*\[[^\[\]]+\]\s*$', '', f).strip()
    f = re.sub(r'(倍|°)$', '', f).strip()
    return f


def same(a, b, n=6):
    a, b = sp.sympify(a), sp.sympify(b)
    syms = sorted((a.free_symbols | b.free_symbols), key=str)
    names = {}
    for x in syms:
        names.setdefault(x.name, []).append(x)
    rnd = random.Random(1)
    for _ in range(n):
        sub = {}
        for nm, xs in names.items():
            v = rnd.uniform(0.3, 1.7) if nm != 'theta' else rnd.uniform(0.1, 1.4)
            for x in xs:
                sub[x] = v
        va, vb = complex(a.subs(sub).evalf()), complex(b.subs(sub).evalf())
        if abs(va - vb) > 1e-9 * max(1, abs(va), abs(vb)):
            return False
    return True


def check(pid, derived, k=0, rel=None):
    f = field(pid, k)
    if rel:
        lhs, rhs = f.split(rel)
        f = rhs
    got = parse(f)
    ok = same(got, derived)
    print(('OK ' if ok else 'NG ') + f'({pid}) 入力用の答え {f}  ／ 導出 {sp.simplify(derived)}')
    assert ok, f'({pid}) が一致しません'


def check_tex(pid, *needles):
    tex = next(a['tex'] for a in PROB['answer']['parts'] if a['id'] == pid)
    for nd in needles:
        assert nd in tex, f'({pid}) の答え {tex} に {nd} がありません'
    print(f'OK ({pid}) 答えの文に ' + '，'.join(needles))


# ----------------------------------------------------------------
import numpy as np
from scipy.integrate import solve_ivp
m, q, Bm, v, th, D = S('m q B v theta D')
vperp, vpar = v * sp.sin(th), v * sp.cos(th)
r = m * vperp / (q * Bm)
T = 2 * sp.pi * r / vperp
check('1', r)
check('2', sp.simplify(T))
check('3', vpar * T)
check_tex('4', '変わらない')
t = D / vpar
phi = 2 * sp.pi / T * t
s = 2 * r * sp.Abs(sp.sin(phi / 2))
check('5', s)
Bmin = sp.solve(sp.Eq(D, vpar * 2 * sp.pi * m / (q * Bm)), Bm)[0]
check('6', Bmin)
# 数値積分で (5) を確かめる
pm = dict(m=1.0, q=1.0, B=2.0, v=1.3, th=0.6, D=2.7)
def f(_, y):
    vel = y[3:]
    a = pm['q'] / pm['m'] * np.cross(vel, [0, 0, pm['B']])
    return np.concatenate([vel, a])
y0 = [0, 0, 0, pm['v'] * np.sin(pm['th']), 0, pm['v'] * np.cos(pm['th'])]
tend = pm['D'] / (pm['v'] * np.cos(pm['th']))
sol = solve_ivp(f, [0, tend], y0, rtol=1e-11, atol=1e-12)
dist = np.hypot(sol.y[0, -1], sol.y[1, -1])
pred = float(s.subs({m: pm['m'], q: pm['q'], Bm: pm['B'], v: pm['v'], th: pm['th'], D: pm['D']}))
assert abs(dist - pred) < 1e-7, (dist, pred)
print('OK 数値積分：スクリーン上の距離', dist)
