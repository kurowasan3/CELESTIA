# 検算：MEC-FRC-105（なめらかな円環に通した小環とばねのつり合い（力の三角形と三角形 OAP の相似））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-FRC-105.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-FRC-105/problem.json').read_text(encoding='utf-8'))
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

m, g, k, R, L0, l = S('m', 'g', 'k', 'R', 'L_0', 'l')
ph = sp.Symbol('φ', positive=True)            # ∠AOP（A から時計回り）
# 座標：O 原点、A=(0,R)、P=(R sinφ, R cosφ)
Px, Py = R*sp.sin(ph), R*sp.cos(ph)
AP = 2*R*sp.sin(ph/2)
ux, uy = (0 - Px)/AP, (R - Py)/AP           # P→A の単位ベクトル
Fs = k*(AP - L0)
# 接線方向（OP に垂直）の力のつり合い
tx, ty = sp.cos(ph), -sp.sin(ph)
tang = sp.simplify(Fs*(ux*tx + uy*ty) + (-m*g)*ty)
# 数値で φ を解き、l・F・N を確認
import math
for vals in [dict(m=1, g=9.8, k=30, R=1, L0=0.8), dict(m=0.5, g=9.8, k=12, R=0.6, L0=0.3), dict(m=2, g=9.8, k=50, R=1.2, L0=1.5)]:
    sub = {m: vals['m'], g: vals['g'], k: vals['k'], R: vals['R'], L0: vals['L0']}
    f = sp.lambdify(ph, tang.subs(sub))
    roots, prev = [], None
    for i in range(1, 3600):
        x = i*math.pi/3600
        if prev is not None and f(x)*prev < 0:
            roots.append(sp.nsolve(tang.subs(sub), ph, x))
        prev = f(x)
    assert len(roots) == 1, roots
    r = roots[0]
    lnum = float(AP.subs(sub).subs(ph, r))
    assert abs(lnum - float((k*R*L0/(k*R - m*g)).subs(sub))) < 1e-9
    Fnum = float(Fs.subs(sub).subs(ph, r))
    assert abs(Fnum - vals['m']*vals['g']*lnum/vals['R']) < 1e-9
    # 半径方向（外向き正）：N + F·(OP 方向) + mg·(OP 方向) = 0
    rx, ry = sp.sin(r), sp.cos(r)
    Nout = -(Fnum*float((ux*rx + uy*ry).subs(sub).subs(ph, r)) + (-vals['m']*vals['g'])*float(ry))
    assert abs(Nout - vals['m']*vals['g']) < 1e-9          # 大きさ mg、外向き（正）
    assert vals['k']*(2*vals['R'] - vals['L0']) > 2*vals['m']*vals['g']
# 条件を満たさない例では最下点以外の解がない
sub = {m: 1, g: 9.8, k: 10, R: 1, L0: 0.5}           # k(2R-L0)=15 < 2mg=19.6
f = sp.lambdify(ph, tang.subs(sub))
assert all(f(i*math.pi/3600)*f((i+1)*math.pi/3600) > 0 for i in range(1, 3599))
check('1', m*g*l/R)
check('2', m*g)
check('3', k*R*L0/(k*R - m*g))
print('(4) OK: k(2R-L0) > 2mg（数値で確認）')
print('すべて一致')
