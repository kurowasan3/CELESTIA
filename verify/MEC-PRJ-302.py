# 検算：MEC-PRJ-302（崖の上から投げた小球の到達距離を最大にする角度（着地速度との直交と，高さ 0 の極限））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-PRJ-302.py
import json, math, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-PRJ-302/problem.json').read_text(encoding='utf-8'))
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

v0, g, h, th = S('v_0', 'g', 'h', 'θ')
t = sp.Symbol('t')
y = v0 * sp.sin(th) * t - g * t**2 / 2
sols = sp.solve(sp.Eq(y, -h), t)
t1 = [s for s in sols if float(s.subs({v0: 1, g: 1, h: 1, th: 0.5})) > 0][0]
check('1', t1, cond=lambda v: v[th] < math.pi / 2)
Vx, Vy = v0 * sp.cos(th), v0 * sp.sin(th) - g * t1
check('2', sp.sqrt(Vx**2 + Vy**2), cond=lambda v: v[th] < math.pi / 2)
R = v0 * sp.cos(th) * t1
# (3) gR = v0 V sinψ（外積で確かめる）
Vv = sp.sqrt(v0**2 + 2 * g * h)
cross = sp.Abs(v0 * sp.cos(th) * Vy - v0 * sp.sin(th) * Vx)
for vals in ({v0: 1.3, g: 2.1, h: 0.7, th: 0.4}, {v0: 2.0, g: 9.8, h: 3.0, th: 1.1}):
    assert abs(float((g * R - cross).subs(vals))) < 1e-9
# (4) 数値最大化
for vals in ({v0: 1.0, g: 1.0, h: 1.0}, {v0: 14.0, g: 9.8, h: 30.0}, {v0: 3.0, g: 9.8, h: 0.2}):
    Rf = sp.lambdify(th, R.subs(vals))
    best = max((Rf(k * (math.pi / 2) / 200000), k) for k in range(1, 200000))
    thb = best[1] * (math.pi / 2) / 200000
    Rmax_expr = grading('4', 0).subs({Sym('v_0'): vals[v0], Sym('g'): vals[g], Sym('h'): vals[h]})
    tan_expr = grading('4', 1).subs({Sym('v_0'): vals[v0], Sym('g'): vals[g], Sym('h'): vals[h]})
    assert abs(best[0] - float(Rmax_expr)) < 1e-6 * best[0], (best[0], Rmax_expr)
    assert abs(math.tan(thb) - float(tan_expr)) < 1e-3, (math.tan(thb), tan_expr)
print('(4) OK: 数値最大化と一致')
check('4', v0 * Vv / g, 0)
check('4', v0 / Vv, 1)
check('5', sp.Rational(14) * sp.sqrt(14**2 + 2 * sp.Rational(98, 10) * 30) / sp.Rational(98, 10))
# (6) 極限
assert sp.limit(v0 / Vv, h, 0) == 1 and sp.limit(v0 / Vv, h, sp.oo) == 0
print('すべて一致')
