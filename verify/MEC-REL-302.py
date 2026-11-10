# 検算：MEC-REL-302（風の中を往復する飛行機の所要時間（風の向きによる違いと，どんな風でも往復は遅くなる理由））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-REL-302.py
import json, math, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-REL-302/problem.json').read_text(encoding='utf-8'))
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

L, c, w, phi = S('L', 'c', 'w', 'φ')
ok = lambda v: v[w] < v[c]
check('1', 2 * L / c)
check('2', L / (c + w) + L / (c - w), cond=ok)
# (3) 機首角 α：垂直成分を打ち消す
al = sp.asin(w / c)
V = c * sp.cos(al)
check('3', 2 * L / V, 0, cond=ok)
check('3', w / c, 1, cond=ok)
# (4) 一般の向き：ベクトルで解く。行き：c(cosβ, sinβ) + w(cosφ, sinφ) が x 軸方向
def leg(sign):
    # sign=+1 行き（+x へ）、-1 帰り（-x へ）
    beta = sp.Symbol('beta')
    vy = c * sp.sin(beta) + w * sp.sin(phi)
    # 垂直成分 0 → sinβ = -w sinφ / c、x 成分 = sign*sqrt(c^2 - w^2 sin^2φ) + w cosφ
    vx = sign * sp.sqrt(c**2 - w**2 * sp.sin(phi)**2) + w * sp.cos(phi)
    return L / (sign * vx)
Tphi = leg(1) + leg(-1)
check('4', Tphi, cond=lambda v: v[w] < v[c])
# 数値シミュレーション（機首角を数値で探す）でも確かめる
import cmath
cv, wv, ph, Lv = 1.0, 0.6, 0.7, 1.0
def ground_speed(dirn):
    # 垂直成分が 0 になる機首角は sin b = -w sinφ / c の2解。そのうち進みたい向きに進む方
    b1 = math.asin(-wv * math.sin(ph) / cv)
    for b in (b1, math.pi - b1):
        vx = cv * math.cos(b) + wv * math.cos(ph)
        vy = cv * math.sin(b) + wv * math.sin(ph)
        assert abs(vy) < 1e-12
        if vx * dirn > 0:
            return vx
Tnum = Lv / ground_speed(1) + Lv / -ground_speed(-1)
Tf = float(Tphi.subs({L: Lv, c: cv, w: wv, phi: ph}))
assert abs(Tnum - Tf) < 1e-3, (Tnum, Tf)
# (5) 極限と大小
assert sp.simplify(Tphi.subs(phi, 0) - (L / (c + w) + L / (c - w))) == 0
assert sp.simplify(Tphi.subs(phi, sp.pi / 2) - 2 * L / sp.sqrt(c**2 - w**2)) == 0
for wv in (0.1, 0.5, 0.9):
    T0, T1, T2 = 2 / 1, 2 / (1 - wv**2), 2 / math.sqrt(1 - wv**2)
    assert T1 > T2 > T0
print('すべて一致')
