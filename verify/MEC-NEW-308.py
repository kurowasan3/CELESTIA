# 検算：MEC-NEW-308（あらい床上に重ねた2物体の上の物体を引く（下の物体が動くかどうかの条件と，加速度の上限））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-NEW-308.py
import json, math, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-NEW-308/problem.json').read_text(encoding='utf-8'))
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

m, M, g, F, m1, m2 = S('m', 'M', 'g', 'F', 'μ_1', 'μ_2')
check('1', F)
check('2', m1 * m * g)
check('3', m2 * (M + m) * g, 0)
a = F / (M + m) - m2 * g
check('3', a, 1)
f = F - m * a
F3 = sp.solve(sp.Eq(f, m1 * m * g), F)[0]
check('4', F3)
check('5', F / m - m1 * g, 0)
check('5', (m1 * m * g - m2 * (M + m) * g) / M, 1)
# 数値シミュレーション：F をゆっくり増やし、各 F で静止・一体・滑りを判定して加速度を比べる
def regime(Fv, mv, Mv, u1, u2, gv=9.8):
    # 1) 全体静止？
    if Fv <= min(u1 * mv * gv, u2 * (mv + Mv) * gv):
        return ('rest', 0.0, 0.0)
    # 2) 一体で動く？
    aa = (Fv - u2 * (mv + Mv) * gv) / (mv + Mv)
    if aa > 0 and Fv - mv * aa <= u1 * mv * gv:
        return ('together', aa, aa)
    # 3) A が滑る
    aA = Fv / mv - u1 * gv
    aB = max(0.0, (u1 * mv * gv - u2 * (mv + Mv) * gv) / Mv)
    return ('slip', aA, aB)
mv, Mv, u1, u2 = 1.0, 2.0, 0.7, 0.2      # μ1 m > μ2 (M+m)
F2v = u2 * (mv + Mv) * 9.8; F3v = float(F3.subs({m: mv, M: Mv, m1: u1, m2: u2, g: 9.8}))
assert regime(F2v * 0.99, mv, Mv, u1, u2)[0] == 'rest'
assert regime((F2v + F3v) / 2, mv, Mv, u1, u2)[0] == 'together'
r = regime(F3v * 1.5, mv, Mv, u1, u2)
assert r[0] == 'slip' and abs(r[2] - (u1 * mv - u2 * (mv + Mv)) * 9.8 / Mv) < 1e-12
assert abs(regime(F3v * 0.999999, mv, Mv, u1, u2)[2] - r[2]) < 1e-4   # F3 で加速度がつながる
mv, Mv, u1, u2 = 1.0, 2.0, 0.3, 0.2      # μ1 m < μ2 (M+m)：B は動かない
for Fv in (1, 5, 50, 500):
    assert regime(Fv, mv, Mv, u1, u2)[2] == 0.0
print('すべて一致')
