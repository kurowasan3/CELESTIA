# 検算：MEC-ENE-104（机の端から滑り落ちる一様な鎖（重心で考える位置エネルギーと，変化する摩擦力の仕事））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-ENE-104.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-ENE-104/problem.json').read_text(encoding='utf-8'))
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

m, g, L, a, x, mu, muk = S('m', 'g', 'L', 'a', 'x', 'μ', 'μk')
U = -(m/L)*x*g*x/2
v1 = sp.sqrt(2*(U.subs(x, a) - U.subs(x, L))/m)
check('1', v1, cond=lambda v: v[a] < v[L])
# (2) 垂れた部分の重力 > 最大静止摩擦力
ac = sp.solve(sp.Eq((m/L)*a*g, mu*(m/L)*(L - a)*g), a)[0]
assert sp.simplify(ac - mu*L/(1 + mu)) == 0
print('(2) OK: a >', ac)
W = -sp.integrate(muk*(m/L)*(L - x)*g, (x, a, L))
assert sp.simplify(W + muk*m*g*(L - a)**2/(2*L)) == 0
print('(3) OK: W =', sp.factor(W))
v2 = sp.sqrt(2*(U.subs(x, a) - U.subs(x, L) + W)/m)
assert same(v2, sp.sqrt(g*((L**2 - a**2) - muk*(L - a)**2)/L), cond=lambda v: v[a] < v[L])
# 運動方程式 m dv/dt = (m/L) x g - μ'(m/L)(L-x) g を数値積分して確認
import math
for (Lv, av, mk) in [(1.0, 0.4, 0.3), (2.0, 0.5, 0.1), (1.5, 1.0, 0.0)]:
    xx, vv, dt = av, 0.0, 1e-6
    while xx < Lv:
        acc = 9.8*(xx - mk*(Lv - xx))/Lv
        vv += acc*dt; xx += vv*dt
    ref = float(v2.subs({g: 9.8, L: Lv, a: av, muk: mk}))
    assert abs(vv - ref) < 1e-3, (vv, ref)
print('(4) OK: 数値積分と一致')
print('すべて一致')
