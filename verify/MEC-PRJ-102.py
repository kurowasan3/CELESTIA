# 検算：MEC-PRJ-102（速度が初速度と垂直になる瞬間（速度ベクトルの三角形で見る斜方投射））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-PRJ-102.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-PRJ-102/problem.json').read_text(encoding='utf-8'))
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

v0, g, th, t = S('v_0', 'g', 'θ', 't')
vx, vy = v0*sp.cos(th), v0*sp.sin(th) - g*t
ok = lambda v: v[th] < 1.57
t1 = sp.solve(sp.Eq(vx*v0*sp.cos(th) + vy*v0*sp.sin(th), 0), t)[0]
check('1', t1, cond=ok)
check('2', sp.sqrt(vx**2 + vy.subs(t, t1)**2), cond=ok)
y1 = (v0*sp.sin(th)*t - g*t**2/2).subs(t, t1)
check('3', y1, cond=ok)
# (4) t1 <= 滞空時間 2 v0 sinθ / g  <=>  θ >= 45°
for d in range(1, 90):
    r = sp.rad(d)
    before = (t1 - 2*v0*sp.sin(th)/g).subs({th: r, v0: 1, g: 1}) <= 1e-15
    assert before == (d >= 45), d
    assert (y1.subs({th: r, v0: 1, g: 1}) >= -1e-15) == (d >= 45)
print('(4) OK: 45° <= θ < 90°')
print('すべて一致')
