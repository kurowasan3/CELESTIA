# 検算：MEC-FRC-104（水平な力で斜面上に物体を止めておける力の範囲（摩擦角と，いくら押しても滑り上がらない条件））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-FRC-104.py
import json, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-FRC-104/problem.json').read_text(encoding='utf-8'))
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

m, g, th, mu, F = S('m', 'g', 'θ', 'μ', 'F')
N = m*g*sp.cos(th) + F*sp.sin(th)
check('1', N)
down = mu < sp.tan(th)
c1 = lambda v: v[th] < 1.57 and v[mu] < sp.tan(v[th])
c2 = lambda v: c1(v) and v[mu]*sp.tan(v[th]) < 1
F1 = sp.solve(sp.Eq(F*sp.cos(th) + mu*N, m*g*sp.sin(th)), F)[0]
F2 = sp.solve(sp.Eq(F*sp.cos(th), m*g*sp.sin(th) + mu*N), F)[0]
check('2', F1, cond=c1)
check('3', F2, cond=c2)
# (4) 摩擦角による形
lam = sp.atan(mu)
assert same(F1, m*g*sp.tan(th - lam), c1) and same(F2, m*g*sp.tan(th + lam), c2)
# 数値シミュレーション的確認：F をふって、必要な静止摩擦力 f_req（上向き正）が |f_req| <= μN となる範囲
import math
for (t0, u0) in [(0.5, 0.3), (0.9, 0.5), (1.2, 0.8)]:
    ok = []
    for k in range(1, 40000):
        Fv = k*0.001
        Nv = math.cos(t0) + Fv*math.sin(t0)
        freq = math.sin(t0) - Fv*math.cos(t0)
        ok.append(abs(freq) <= u0*Nv)
    lo = (ok.index(True) + 1)*0.001
    hi = None if ok[-1] else (len(ok) - ok[::-1].index(True))*0.001
    assert abs(lo - float(F1.subs({m: 1, g: 1, th: t0, mu: u0}))) < 2e-3
    if u0*math.tan(t0) < 1:
        assert abs(hi - float(F2.subs({m: 1, g: 1, th: t0, mu: u0}))) < 2e-3
    else:
        assert hi is None
print('(4) OK: 範囲を数値で確認、μtanθ>=1 で上限なし')
print('すべて一致')
