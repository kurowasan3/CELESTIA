# 検算：MEC-NEW-305（あらい斜面上の物体と滑車でつるしたおもり（静止する範囲，滑り上がる・滑り下りるときの加速度））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-NEW-305.py
import json, math, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-NEW-305/problem.json').read_text(encoding='utf-8'))
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

M, m, g, th, mu = S('M', 'm', 'g', 'θ', 'μ')
a, T = sp.symbols('a T')
up = sp.solve([sp.Eq(M * a, T - M * g * sp.sin(th) - mu * M * g * sp.cos(th)), sp.Eq(m * a, m * g - T)], [a, T])
check('2', up[a])
vals = {M: 2, m: 3, g: sp.Rational(98, 10), th: sp.pi / 6, mu: sp.Rational(1, 5)}
av = float(up[a].subs(vals)); Tv = float(up[T].subs(vals))
check('3', sig(av), 0); check('3', sig(Tv), 1)
dn = sp.solve([sp.Eq(M * a, M * g * sp.sin(th) - T - mu * M * g * sp.cos(th)), sp.Eq(m * a, T - m * g)], [a, T])
check('4', dn[a])
# 図の力の矢印（T, Mg, N, f）の合力が斜面上向きで大きさ M a になっていることを確認
N = 19.6 * math.cos(math.pi / 6); f = 0.2 * N
Fx = Tv * math.cos(math.pi / 6) - f * math.cos(math.pi / 6) + N * math.cos(math.radians(120))
Fy = Tv * math.sin(math.pi / 6) - f * math.sin(math.pi / 6) + N * math.sin(math.radians(120)) - 19.6
assert abs(math.hypot(Fx, Fy) - 2 * av) < 1e-9 and abs(math.degrees(math.atan2(Fy, Fx)) - 30) < 1e-6
print('すべて一致')
