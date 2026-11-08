# 検算スクリプト（problem.json の入力用の答えを sympy で独立に導いた値と照合する）
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ESY-ROD-101'
PROB = json.loads((Path(__file__).resolve().parent.parent / 'content' / 'problems' / PID[:3] / PID / 'problem.json').read_text(encoding='utf-8'))
GRADING = {a['id']: a.get('grading') for a in PROB['answer']['parts']}
GREEK = {'μ': 'mu', 'θ': 'theta', 'ω': 'omega', 'Φ': 'Phi', 'τ': 'tau', 'Δ': 'Delta', 'φ': 'phi'}
FUNCS = {'sqrt', 'sin', 'cos', 'tan', 'log', 'exp', 'pi'}
TR = standard_transformations + (implicit_multiplication, convert_xor)


def S(names):
    """記号（名前は入力用の答えと同じ：μ_0 → mu_0）"""
    return sp.symbols(names, positive=True)


def parse(src):
    s = src.strip()
    for k, v in GREEK.items():
        s = s.replace(k, ' ' + v + ' ' if v != 'mu' else v)
    s = s.replace('π', ' pi ')
    s = re.sub(r'\s+_', '_', s)
    s = re.sub(r'(\d)\s*([A-Za-z(])', r'\1*\2', s)
    names = set(re.findall(r'[A-Za-z][A-Za-z0-9_]*', s)) - FUNCS
    loc = {n: sp.Symbol(n, positive=True) for n in names}
    loc.update({'sqrt': sp.sqrt, 'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'log': sp.log, 'exp': sp.exp, 'pi': sp.pi})
    return parse_expr(s, local_dict=loc, transformations=TR)


def field(pid, k=0):
    g = GRADING[pid]
    assert g, f'({pid}) に入力用の答えがありません'
    f = [x.strip() for x in g.split(' ; ')][k]
    if '：' in f:
        f = f.split('：', 1)[1]
    f = re.sub(r'\s*\[[^\[\]]+\]\s*$', '', f).strip()
    f = re.sub(r'倍$', '', f).strip()
    return f


def same(a, b, n=6):
    a, b = sp.sympify(a), sp.sympify(b)
    syms = sorted((a.free_symbols | b.free_symbols), key=str)
    names = {}
    for x in syms:
        names.setdefault(x.name, []).append(x)
    rnd = random.Random(1)
    for _ in range(n):
        sub = {}
        for nm, xs in names.items():
            v = rnd.uniform(0.3, 1.7) if nm != 'theta' else rnd.uniform(0.1, 1.4)
            for x in xs:
                sub[x] = v
        va, vb = complex(a.subs(sub).evalf()), complex(b.subs(sub).evalf())
        if abs(va - vb) > 1e-9 * max(1, abs(va), abs(vb)):
            return False
    return True


def check(pid, derived, k=0, rel=None):
    f = field(pid, k)
    if rel:
        lhs, rhs = f.split(rel)
        f = rhs
    got = parse(f)
    ok = same(got, derived)
    print(('OK ' if ok else 'NG ') + f'({pid}) 入力用の答え {f}  ／ 導出 {sp.simplify(derived)}')
    assert ok, f'({pid}) が一致しません'


def check_tex(pid, *needles):
    tex = next(a['tex'] for a in PROB['answer']['parts'] if a['id'] == pid)
    for nd in needles:
        assert nd in tex, f'({pid}) の答え {tex} に {nd} がありません'
    print(f'OK ({pid}) 答えの文に ' + '，'.join(needles))

# ----------------------------------------------------------------
import numpy as np
from scipy.integrate import solve_ivp
B, l, v, R, m, g, th = S('B l v R m g theta')
# 座標：x 水平（斜面を下る向き）、y 鉛直上向き、z は棒の向き。下る向き d = (cosθ, -sinθ)、外向き法線 n = (sinθ, cosθ)
vel = v * sp.Matrix([sp.cos(th), -sp.sin(th), 0])
Bv = sp.Matrix([0, B, 0])
emf_dir = vel.cross(Bv)            # 単位電荷の受ける力（棒の向き z に沿う成分が起電力/長さ）
V = sp.simplify(sp.Abs(emf_dir[2]) * l)
check('1', V)
I = V / R
# 電流は起電力の向き：z 成分の符号に合わせる
Ivec = sp.sign(emf_dir[2]) * I * sp.Matrix([0, 0, 1])
Fm = (Ivec * l).cross(Bv)
Fm = sp.simplify(Fm.subs(sp.sign(emf_dir[2]), -1)) if emf_dir[2].has(sp.sign) else Fm
Fm = sp.simplify(Fm)
assert sp.simplify(Fm[1]) == 0
assert sp.simplify(Fm[0] + B**2 * l**2 * v * sp.cos(th) / R) == 0   # 水平で -x（斜面を上る側）
check('2', sp.Abs(Fm[0]))
dvec = sp.Matrix([sp.cos(th), -sp.sin(th), 0]); nvec = sp.Matrix([sp.sin(th), sp.cos(th), 0])
grav = sp.Matrix([0, -m * g, 0])
a = sp.simplify((grav + Fm).dot(dvec) / m)
check('3', a)
vf = sp.solve(sp.Eq(a, 0), v)[0]
check('4', vf)
N = sp.simplify(-(grav + Fm).dot(nvec).subs(v, vf))
check('5', N)
P = sp.simplify((I**2 * R).subs(v, vf))
check('6', P)
assert sp.simplify(P - m * g * sp.sin(th) * vf) == 0
# 数値：終端速度に近づく
p = {B: 0.8, l: 0.5, R: 1.5, m: 0.2, g: 9.8, th: 0.5}
af = sp.lambdify(v, a.subs(p))
s = solve_ivp(lambda t, y: [af(y[0])], [0, 120], [0.0], rtol=1e-11, atol=1e-12)
assert abs(s.y[0][-1] - float(vf.subs(p))) < 1e-6
print('OK 数値：終端速度', s.y[0][-1])
# 図の検査：解答図の力のつり合いの比 N : mg = 1/cosθ（θ = 30°）
import json as _j
spec = _j.loads((Path(__file__).resolve().parent.parent / 'content/problems/ESY/ESY-ROD-101/a-forces.fig.json').read_text())
arr = {i['phys']: i for i in spec['items'] if i.get('t') == 'arrow' and i.get('phys')}
assert abs(arr['force:N']['len'] / arr['force:mg']['len'] - 1 / np.cos(np.radians(30))) < 1e-6
assert abs(arr['force:F']['len'] / arr['force:mg']['len'] - np.tan(np.radians(30))) < 1e-6
print('OK 解答図の矢印の長さの比')
