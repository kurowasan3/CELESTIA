# 検算：MEC-KIN-302（記録タイマーのテープから平均の速度・加速度・最初の打点での速さを求める）
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-KIN-302.py
import json, math, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-KIN-302/problem.json').read_text(encoding='utf-8'))
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

L = [sp.Rational(x) / 100 for x in ('1.5', '3.5', '5.5', '7.5', '9.5')]   # m
T = sp.Rational(5, 50)
check('1', T)
vbar = [l / T for l in L]
check('2', vbar[2])
d = [L[i + 1] - L[i] for i in range(4)]
assert len(set(d)) == 1
a = d[0] / T**2
check('3', a)
# 等加速度の式で全区間を再現して、最初の打点の速さ v0 を最小二乗なしで決める
v0, aa = sp.symbols('v0 aa')
sol = sp.solve([sp.Eq(v0 * T + aa * T**2 / 2, L[0]), sp.Eq(v0 * 2 * T + aa * (2 * T)**2 / 2, L[0] + L[1])], [v0, aa])
assert sol[aa] == a
for k in range(1, 6):
    assert sol[v0] * k * T + a * (k * T)**2 / 2 == sum(L[:k])
check('4', sol[v0] + a * 2 * T)
check('5', sol[v0])
print('すべて一致')
