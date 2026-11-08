# 検算：MEC-REL-101（落下する小球と投げ上げた小球の空中での出会い（相対運動で見る））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-REL-101.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-REL-101/problem.json').read_text(encoding='utf-8'))
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

H, v0, g, t = S('H', 'v_0', 'g', 't')
yA = H - g*t**2/2
yB = v0*t - g*t**2/2
assert sp.diff(yB - yA, t, 2) == 0          # 相対加速度 0
check('1', sp.diff(yB, t, 2) - sp.diff(yA, t, 2))
t1 = sp.solve(sp.Eq(yA, yB), t)[0]
check('2', t1)
# (3) 上昇中：B の速度 v0 - g t1 > 0  <=>  v0^2 > g H
cond = sp.simplify(v0 - g*t1)
assert sp.simplify(cond - (v0**2 - g*H)/v0) == 0
num = {H: sp.Rational(196, 10), v0: sp.Rational(196, 10), g: sp.Rational(98, 10)}
check('4', yA.subs(t, t1).subs(num))
assert (v0**2 - g*H).subs(num) > 0
print('すべて一致')
