# 検算：MEC-NEW-307（坂を下る電車の天井からつるしたおもりの傾き（加速度計として使う・糸が床に垂直になる意味））
# problem.json の答え（grading）を、独立に導いた結果と照合する。python3 -I verify/MEC-NEW-307.py
import json, math, random, re, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

P = json.loads((pathlib.Path(__file__).resolve().parent.parent / 'content/problems/MEC/MEC-NEW-307/problem.json').read_text(encoding='utf-8'))
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

m, g, al, a = sp.symbols('m g α a', positive=True)
T, psi = sp.symbols('T psi')
# 坂に沿って下向き x、床に垂直 y。張力の成分 (-T sinψ, T cosψ)
Tx, Ty = sp.symbols('Tx Ty')
sol = sp.solve([sp.Eq(m * a, m * g * sp.sin(al) + Tx), sp.Eq(0, Ty - m * g * sp.cos(al))], [Tx, Ty])
tanpsi = -sol[Tx] / sol[Ty]
Tmag = sp.sqrt(sol[Tx]**2 + sol[Ty]**2)
check('1', Tmag.subs(a, 0))
check('2', tanpsi, 0)
check('2', Tmag, 1, cond=lambda v: True)
gg, sa, ca = 9.8, 0.10, 0.995
check('3', sig(gg * sa - gg * ca * 0.050))
check('4', sig((gg * sa + 0.98) / (gg * ca)))
aperp = sp.solve(sp.Eq(tanpsi, 0), a)[0]
assert sp.simplify(aperp - g * sp.sin(al)) == 0
check('5', sp.simplify(Tmag.subs(a, aperp)), cond=lambda v: v[al] < math.pi / 2)
# 図（α=30°, a=g sinα/2）：張力の向き 76.1°
Tv = (sol[Tx].subs({m: 1, g: 1, al: sp.pi / 6, a: sp.Rational(1, 4)}), sol[Ty].subs({m: 1, g: 1, al: sp.pi / 6, a: sp.Rational(1, 4)}))
# 坂座標 → 水平・鉛直座標（x 下り方向 = (cos30, -sin30), y 法線 = (sin30, cos30)）
hx = float(Tv[0]) * math.cos(math.pi / 6) + float(Tv[1]) * math.sin(math.pi / 6)
hy = -float(Tv[0]) * math.sin(math.pi / 6) + float(Tv[1]) * math.cos(math.pi / 6)
assert abs(math.degrees(math.atan2(hy, hx)) - 76.1) < 0.1, math.degrees(math.atan2(hy, hx))
print('すべて一致')
