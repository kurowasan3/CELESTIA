# 検算：MEC-FRC-307（2つの滑車にかけた糸の結び目につるしたおもり（力の三角形とつり合いが存在する条件））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-FRC-307.py
import json, math, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-FRC-307/problem.json').read_text(encoding='utf-8'))
ANS = {a['id']: a for a in P['answer']['parts']}
FUN = {'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'sqrt': sp.sqrt, 'log': sp.log, 'π': sp.pi, 'pi': sp.pi}

def grading(pid, k=0):
    s = ANS[pid]['grading'].split(' ; ')[k].strip()
    s = re.sub(r'^[^：]{1,80}：\s*', '', s)
    s = re.sub(r'\s*\[[^\]]+\]$', '', s).replace('倍', '').replace('°', '').replace('^', '**')
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
    for _ in range(4000):
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

def sig(x, n=2):
    """有効数字 n 桁に丸めた値（答えに書いた値と同じもの）"""
    x = float(x)
    if x == 0: return sp.Integer(0)
    e = math.floor(math.log10(abs(x))) - n + 1
    return sp.Rational(round(x / 10**e)) * sp.Rational(10)**e

S = lambda *names: sp.symbols(names, positive=True)
Sym = lambda name: sp.Symbol(name, positive=True)

m1, m2, m3 = S('m_1', 'm_2', 'm_3')
a_, b_ = sp.symbols('a b')
ca = (m3**2 + m1**2 - m2**2) / (2 * m1 * m3)
cb = (m3**2 + m2**2 - m1**2) / (2 * m2 * m3)
# つり合いの2式を満たすか（三角形が成立する範囲で）
okc = lambda v: abs(v[m1] - v[m2]) < v[m3] < v[m1] + v[m2]
assert same(m1 * sp.sqrt(1 - ca**2) - m2 * sp.sqrt(1 - cb**2), 0, okc)
assert same(m1 * ca + m2 * cb - m3, 0, okc)
check('2', ca, 0, cond=okc); check('2', cb, 1, cond=okc)
# (3) 数値
c1 = ca.subs({m1: 3, m2: 4, m3: 5}); c2 = cb.subs({m1: 3, m2: 4, m3: 5})
ang = sp.acos(c1) + sp.acos(c2)
assert abs(float(ang) - math.pi / 2) < 1e-12 and c1 * c2 - sp.sqrt(1 - c1**2) * sp.sqrt(1 - c2**2) == 0
check('3', 90, 0)
d, y = S('d', 'y')
ta, tb = sp.sqrt(1 - c1**2) / c1, sp.sqrt(1 - c2**2) / c2
check('3', sp.solve(sp.Eq(y * ta + y * tb, d), y)[0], 1)
# (5) 存在条件：ランダムな質量で、cos が (0,1) に入る ⇔ sqrt|m1^2-m2^2| < m3 < m1+m2
rnd = random.Random(3)
for _ in range(20000):
    a, b, c = (rnd.uniform(0.1, 3) for _ in range(3))
    x = (c * c + a * a - b * b) / (2 * a * c); z = (c * c + b * b - a * a) / (2 * b * c)
    exists = 0 < x < 1 and 0 < z < 1
    cond = math.sqrt(abs(a * a - b * b)) < c < a + b
    assert exists == cond, (a, b, c, x, z)
print('(5) OK: 存在条件')
print('すべて一致')
