# 検算：MEC-NEW-306（水平に引く力の F-t グラフから，あらい床の上の物体の v-t グラフをかき，止まる時刻と移動距離を求める）
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-NEW-306.py
import json, math, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-NEW-306/problem.json').read_text(encoding='utf-8'))
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

m, g = 2.0, 9.8
fs, fk = 0.50 * m * g, 0.25 * m * g
def F(t): return 14.7 if t < 2.0 else (4.9 if t < 4.0 else 0.0)
assert F(0) > fs
# 数値で運動を追う
dt = 1e-5; t = x = v = 0.0; tstop = None; v2 = None
while True:
    if v > 0:
        a = (F(t) - fk) / m
    else:
        a = (F(t) - fk) / m if F(t) > fs else 0.0
    v_new = v + a * dt
    if v > 0 and v_new <= 0 and t > 4:
        tstop = t; break
    x += (v + v_new) / 2 * dt; v = v_new; t += dt
    if v2 is None and t >= 2.0: v2 = v
check('1', sp.Rational(147, 10) / 2 - sp.Rational(49, 10) / 2)
assert abs(v2 - 9.8) < 1e-3
check('2', sp.Rational(98, 10))
assert abs(tstop - 8.0) < 1e-3, tstop
check('3', 8)
assert abs(x - 49) < 1e-3, x
check('4', 49)
check('5', 0)
print('すべて一致')
