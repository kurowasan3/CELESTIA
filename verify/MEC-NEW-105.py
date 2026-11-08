# 検算：MEC-NEW-105（動滑車でつながれた机上の物体とおもり（糸の長さ一定から決まる加速度の関係））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-NEW-105.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-NEW-105/problem.json').read_text(encoding='utf-8'))
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

m, M, g, T, aB, mu, muk = S('m', 'M', 'g', 'T', 'a_B', 'μ', 'μk')
x = sp.Symbol('x', positive=True)
# 拘束：糸の全長 = (机上の部分) + (定滑車から動滑車まで) + (動滑車から C まで)
# A が x 動くと机上の部分が x 減る → 下の2本が合わせて x 伸びる → 動滑車は x/2 下がる
y = sp.solve(sp.Eq(-x + 2*sp.Symbol('y'), 0), sp.Symbol('y'))[0]
check('1', x/y)
sol = sp.solve([sp.Eq(m*2*aB, T), sp.Eq(M*aB, M*g - 2*T)], [aB, T], dict=True)[0]
check('2', sol[aB])
check('3', sol[T])
# エネルギーで確認
h = sp.Symbol('h', positive=True); vB = sp.sqrt(2*sol[aB]*h)
assert sp.simplify(M*g*h - (M*vB**2/2 + m*(2*vB)**2/2)) == 0
# (4)
Tstatic = M*g/2
assert sp.simplify(sp.solve(sp.Eq(Tstatic, mu*m*g), M)[0] - 2*mu*m) == 0
sol2 = sp.solve([sp.Eq(m*2*aB, T - muk*m*g), sp.Eq(M*aB, M*g - 2*T)], [aB, T], dict=True)[0]
assert sp.simplify(sol2[aB] - (M - 2*muk*m)*g/(M + 4*m)) == 0
print('(4) OK: M > 2μm, a_B =', sp.factor(sol2[aB]))
print('すべて一致')
