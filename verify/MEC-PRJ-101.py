# 検算：MEC-PRJ-101（斜面に向けた斜方投射（斜面上の到達距離と，それを最大にする投射角））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-PRJ-101.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-PRJ-101/problem.json').read_text(encoding='utf-8'))
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

v0, g, th, al, t = S('v_0', 'g', 'θ', 'α', 't')
# 水平・鉛直の座標で独立に解く
x = v0*sp.cos(th)*t; y = v0*sp.sin(th)*t - g*t**2/2
Tsol = [s for s in sp.solve(sp.Eq(y, x*sp.tan(al)), t) if s != 0][0]
ok = lambda v: 0 < v[al] < 1.4 and v[al] < v[th] < sp.pi/2
check('1', Tsol, cond=ok)
L = x.subs(t, Tsol)/sp.cos(al)
check('2', L, cond=ok)
# (3) dL/dθ = 0
dL = sp.diff(L, th)
for a0 in [0.2, 0.5, 0.9]:
    f = sp.lambdify(th, L.subs({al: a0, v0: 1, g: 1}))
    grid = max((f(a0 + k*(1.5707963 - a0)/20000), a0 + k*(1.5707963 - a0)/20000) for k in range(1, 20000))
    assert abs(grid[1] - (float(sp.pi)/4 + a0/2)) < 1e-3
    thopt = sp.nsolve(dL.subs({al: a0, v0: 1, g: 1}), th, float(sp.pi/4 + a0/2) + 0.05)
    assert abs(thopt - (sp.pi/4 + a0/2)) < 1e-9
    Lmax = L.subs({al: a0, v0: 1, g: 1, th: thopt})
    assert abs(Lmax - 1/(1 + sp.sin(a0))) < 1e-9
check('3', sp.pi/4 + al/2, 0)
check('3', L.subs(th, sp.pi/4 + al/2), 1, cond=lambda v: 0 < v[al] < 1.4)
# (4) 着地時の速度が斜面に垂直：速度と斜面方向の内積 0
vx = v0*sp.cos(th); vy = v0*sp.sin(th) - g*Tsol
dot = sp.simplify(vx*sp.cos(al) + vy*sp.sin(al))
for a0 in [0.2, 0.5, 0.9]:
    th0 = a0 + sp.atan(1/(2*sp.tan(a0)))
    assert abs(dot.subs({al: a0, th: th0, v0: 1.3, g: 2.1}).evalf()) < 1e-12
print('(4) OK: tan(θ-α)=1/(2tanα)')
print('すべて一致')
