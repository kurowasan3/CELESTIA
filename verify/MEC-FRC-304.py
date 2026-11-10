# 検算：MEC-FRC-304（おもりの質量とばねの長さの測定グラフからばね定数と自然の長さを読む）
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-FRC-304.py
import json, math, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-FRC-304/problem.json').read_text(encoding='utf-8'))
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

import numpy as np
m = np.array([0, 50, 100, 150, 200, 250]) / 1000      # kg
L = np.array([12.0, 14.6, 16.9, 19.6, 21.9, 24.5]) / 100  # m
slope, icpt = np.polyfit(m, L, 1)
# 図の直線は (0, 12.0 cm) と (250 g, 24.5 cm) を通る。最小二乗の直線とほぼ一致することを確認
assert abs(icpt - 0.120) < 0.002 and abs(slope - 0.125 / 0.25) < 0.01, (icpt, slope)
L0 = sp.Rational(12)
check('1', L0)
k = sp.Rational(25, 100) * sp.Rational(98, 10) / (sp.Rational(245, 1000) - sp.Rational(120, 1000))
check('2', sig(k))
mm = (sp.Rational(20) - L0) / (sp.Rational(125, 10) / 250)    # g
check('3', mm)
check('5', sig(2 * k))
print('すべて一致')
