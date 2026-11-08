# 検算：MEC-ENE-103（ばねで打ち出した物体とあらい区間の往復（摩擦による損失は通った道のりに比例））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-ENE-103.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-ENE-103/problem.json').read_text(encoding='utf-8'))
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

k, m, g, mu, d, x0 = 98, sp.Rational(1, 2), sp.Rational(98, 10), sp.Rational(1, 4), sp.Rational(3, 10), sp.Rational(1, 10)
E0 = k*x0**2/2
fr = mu*m*g
check('1', sp.sqrt(2*E0/m))
E1 = E0 - fr*d
check('2', sp.sqrt(2*E1/m))
check('3', sp.sqrt(2*E1/k))
s = E0/fr                              # あらい区間での道のりの合計
assert s == sp.Rational(2, 5)
# 道のりを往復でたどる：右向きに d 通過 → 戻って (s - d) 進む（< d なら区間内で停止）
rem = s - d
assert 0 < rem < d
check('4', rem)
# 数値シミュレーションで確認（時間刻み）
import math
P, Q = 0.5, 0.8                        # 位置（ばね1の自然長の端を 0）
x, v, dt, t = -0.1, 0.0, 1e-6, 0
L2 = 1.0                               # ばね2の自然長の端
while t < 20:
    F = 0.0
    if x < 0: F += -98*x
    if x > L2: F += -98*(x - L2)
    if P <= x <= Q and v != 0: F += -math.copysign(float(fr), v)
    vn = v + F/0.5*dt
    if P <= x <= Q and v != 0 and vn*v < 0: v = 0.0; break
    v = vn; x += v*dt; t += dt
assert abs((Q - x) - 0.1) < 1e-3, x
print('数値計算で停止位置 Q から', round(Q - x, 4), 'm')
print('すべて一致')
