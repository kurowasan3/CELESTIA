# 検算：MEC-NEW-104（引かれる板と上に載った物体（一緒に動く条件と，滑り出した後に板から落ちるまでの時間））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-NEW-104.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-NEW-104/problem.json').read_text(encoding='utf-8'))
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

M, m, g, F, mu, muk, L, t = S('M', 'm', 'g', 'F', 'μ', 'μk', 'L', 't')
a = F/(M + m)
check('1', m*a)
Fc = sp.solve(sp.Eq(m*a, mu*m*g), F)[0]
assert sp.simplify(Fc - mu*(M + m)*g) == 0
print('(2) OK: F <=', Fc)
aA = sp.solve(sp.Eq(m*sp.Symbol('x'), muk*m*g), sp.Symbol('x'))[0]
aB = sp.solve(sp.Eq(M*sp.Symbol('x'), F - muk*m*g), sp.Symbol('x'))[0]
assert sp.simplify(aA - muk*g) == 0 and sp.simplify(aB - (F - muk*m*g)/M) == 0
# (4) 床から見た変位の差が L になる時刻
t1 = [s for s in sp.solve(sp.Eq(aB*t**2/2 - aA*t**2/2, L), t) if s.could_extract_minus_sign() is False][0]
assert same(t1, sp.sqrt(2*M*L/(F - muk*(M + m)*g)), cond=lambda v: v[F] > v[muk]*(v[M] + v[m])*v[g])
# 数値の運動で確認（オイラー法）
vals = dict(M=2.0, m=1.0, g=9.8, F=30.0, mu=0.4, muk=0.3, L=1.0)
assert vals['F'] > vals['mu']*(vals['M'] + vals['m'])*vals['g']
xa = xb = va = vb = 0.0; dt = 1e-6; tt = 0.0
while xb - xa < vals['L']:
    va += vals['muk']*vals['g']*dt; vb += (vals['F'] - vals['muk']*vals['m']*vals['g'])/vals['M']*dt
    xa += va*dt; xb += vb*dt; tt += dt
ref = float(t1.subs({M: 2, m: 1, g: 9.8, F: 30, muk: 0.3, L: 1}))
assert abs(tt - ref) < 1e-4, (tt, ref)
print('(3)(4) OK: t1 =', ref)
print('すべて一致')
