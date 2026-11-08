# 検算：MEC-FRC-103（斜め上に引く力と垂直抗力・最大摩擦力（動かすのに必要な力を最小にする引く角度））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-FRC-103.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-FRC-103/problem.json').read_text(encoding='utf-8'))
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

m, g, F, th, mu = S('m', 'g', 'F', 'θ', 'μ')
N = sp.solve(sp.Eq(sp.Symbol('N') + F*sp.sin(th), m*g), sp.Symbol('N'))[0]
check('1', N)
F0 = sp.solve(sp.Eq(F*sp.cos(th), mu*N), F)[0]
check('2', F0, cond=lambda v: v[th] < 1.57)
# (3) μ = 3/4 で dF0/dθ = 0
f = F0.subs(mu, sp.Rational(3, 4))
crit = sp.solve(sp.diff(f, th), th)
crit = [c for c in crit if 0 < float(c) < float(sp.pi/2)]
assert len(crit) == 1
check('3', sp.tan(crit[0]), 0)
check('3', f.subs(th, crit[0]), 1)
# 境界より小さいことの確認
assert float(f.subs({th: crit[0], m: 1, g: 1})) < float(f.subs({th: 0, m: 1, g: 1}))
print('すべて一致')
