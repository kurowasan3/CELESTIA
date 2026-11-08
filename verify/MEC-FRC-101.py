# 検算：MEC-FRC-101（ばねの弾性力とフックの法則（直列につないだばねの伸び））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-FRC-101.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-FRC-101/problem.json').read_text(encoding='utf-8'))
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

m, g, x = sp.Rational(1, 2), sp.Rational(98, 10), sp.Rational(49, 1000)
k = m*g/x
check('1', k)
check('2', 1*g/k*100)               # cm
# (3) 直列：各ばねの張力 = mg、伸びの和
check('3', 2*(m*g/k)*100)
# 合成ばね定数 k/2 で見ても同じ
assert m*g/(k/2)*100 == 2*(m*g/k)*100
print('すべて一致')
