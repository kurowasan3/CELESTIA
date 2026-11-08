# 検算：MEC-ENE-105（あらい曲面に沿ってゆっくり引き上げる仕事（摩擦の仕事は水平距離だけで決まる））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-ENE-105.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-ENE-105/problem.json').read_text(encoding='utf-8'))
ANS = {a['id']: a for a in P['answer']['parts']}
FUN = {'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'sqrt': sp.sqrt, 'log': sp.log, 'π': sp.pi, 'pi': sp.pi}

def grading(pid, k=0):
    s = ANS[pid]['grading'].split(' ; ')[k].strip()
    s = re.sub(r'^[^：]{1,80}：\s*', '', s)
    s = re.sub(r'\s*\[[^\]]+\]$', '', s).replace('倍', '').replace('^', '**')
    s = re.sub(r'_\{(\w+)\}', r'_\1', s)
    names = set(re.findall(r'[A-Za-zα-ωΑ-Ω][A-Za-z0-9_α-ωΑ-Ω]*', s)) - set(FUN)
    loc = dict(FUN); loc.update({n: sp.Symbol(n, positive=True) for n in names})
    return parse_expr(s, local_dict=loc, transformations=standard_transformations + (implicit_multiplication, convert_xor))

def same(e, f, cond=lambda v: True, n=8):
    e, f = sp.sympify(e), sp.sympify(f)
    syms = sorted(e.free_symbols | f.free_symbols, key=str)
    if not syms:
        a, b = complex(e.evalf()), complex(f.evalf())
        return abs(a - b) <= 1e-9 * max(1, abs(a))
    rnd, ok = random.Random(1), 0
    for _ in range(2000):
        v = {s: rnd.uniform(0.1, 3) for s in syms}
        if not cond(v): continue
        a, b = complex(e.subs(v).evalf()), complex(f.subs(v).evalf())
        if abs(a.imag) > 1e-12 or abs(b.imag) > 1e-12: continue
        assert abs(a - b) <= 1e-9 * max(1, abs(a)), (e, f, v)
        ok += 1
        if ok >= n: return True
    raise AssertionError(f'比較できる点が足りない {e} {f}')

def check(pid, expected, k=0, cond=lambda v: True):
    got = grading(pid, k)
    assert same(got, expected, cond), f'({pid}) 不一致: problem.json={got}  検算={sp.simplify(expected)}'
    print(f'({pid}) OK: {got}')

S = lambda *names: sp.symbols(names, positive=True)

m, g, h, d, mu, x = S('m', 'g', 'h', 'd', 'μ', 'x')
# 任意の単調な坂 y = h*s(x/d) に沿ってゆっくり引く：F = mg sinφ + μ mg cosφ、ds = dx/cosφ
import math
def work(shape, m_, g_, h_, d_, mu_, n=200000):
    W = 0.0
    for i in range(n):
        x0, x1 = d_*i/n, d_*(i + 1)/n
        dy = h_*(shape(x1/d_) - shape(x0/d_)); dx = x1 - x0
        ds = math.hypot(dx, dy); c, s = dx/ds, dy/ds
        W += (m_*g_*s + mu_*m_*g_*c)*ds
    return W
for shape in [lambda u: u, lambda u: u*u, lambda u: 1 - (1 - u)**2, lambda u: (1 - math.cos(math.pi*u))/2]:
    assert abs(work(shape, 2.0, 9.8, 3.75, 5.0, 0.25) - 98.0) < 1e-6
check('2', m*g*(h + mu*d))
num = {m: 2, g: sp.Rational(98, 10), h: sp.Rational(375, 100), d: 5, mu: sp.Rational(1, 4)}
check('3', (m*g*(h + mu*d)).subs(num))
assert h.subs(num)/d.subs(num) > mu.subs(num)
check('4', sp.sqrt(2*g*(h - mu*d)).subs(num))
print('(1) OK: 坂の形を4通り変えても W は同じ')
print('すべて一致')
