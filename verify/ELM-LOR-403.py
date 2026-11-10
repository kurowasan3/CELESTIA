# 検算スクリプト（problem.json の入力用の答えを sympy で独立に導いた値と照合する）
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ELM-LOR-403'
PROB = json.loads((Path(__file__).resolve().parent.parent / 'content' / 'problems' / PID[:3] / PID / 'problem.json').read_text(encoding='utf-8'))
GRADING = {a['id']: a.get('grading') for a in PROB['answer']['parts']}
GREEK = {'μ': 'mu', 'θ': 'theta', 'ω': 'omega', 'Φ': 'Phi', 'τ': 'tau', 'Δ': 'Delta', 'φ': 'phi'}
FUNCS = {'abs', 'sqrt', 'sin', 'cos', 'tan', 'log', 'exp', 'pi'}
TR = standard_transformations + (implicit_multiplication, convert_xor)


def S(names):
    """記号（名前は入力用の答えと同じ：μ_0 → mu_0）"""
    return sp.symbols(names, positive=True)


def parse(src):
    s = src.strip()
    for k, v in GREEK.items():
        s = s.replace(k, ' ' + v + ' ' if v != 'mu' else v)
    s = s.replace('π', ' pi ').replace('×', '*')
    s = re.sub(r'\s+_', '_', s)
    s = re.sub(r'(\d)\s*([A-Za-z(])', r'\1*\2', s)
    names = set(re.findall(r'[A-Za-z][A-Za-z0-9_]*', s)) - FUNCS
    loc = {n: sp.Symbol(n, positive=True) for n in names}
    loc.update({'abs': sp.Abs, 'sqrt': sp.sqrt, 'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'log': sp.log, 'exp': sp.exp, 'pi': sp.pi})
    return parse_expr(s, local_dict=loc, transformations=TR)


def field(pid, k=0):
    g = GRADING[pid]
    assert g, f'({pid}) に入力用の答えがありません'
    f = [x.strip() for x in g.split(' ; ')][k]
    if '：' in f:
        f = f.split('：', 1)[1]
    f = re.sub(r'\s*\[[^\[\]]+\]\s*$', '', f).strip()
    f = re.sub(r'(倍|°)$', '', f).strip()
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
m, q, Bm, v, a = S('m q B v a')
r = m * v / (q * Bm)
check('1', r)
check('2', a / r)
B0 = sp.solve(sp.Eq(r, a), Bm)[0]
check('3', B0)
# 1/4 周
T0 = 2 * sp.pi * m / (q * B0)
check('4', sp.simplify(T0 / 4))
# (5)(6) 数値積分：B 紙面の裏から表(+z)、q>0
def run(p0, v0, aa=1.0, vv=1.0):
    # m = q = 1, B0 = v/a
    Bz = vv / aa
    def f(_, y):
        return [y[2], y[3], y[3] * Bz, -y[2] * Bz]
    def leave(_, y):
        return y[0]**2 + y[1]**2 - aa**2 * (1 + 1e-9)
    leave.terminal, leave.direction = True, 1
    sol = solve_ivp(f, [0, 20], [*p0, *v0], events=leave, rtol=1e-11, atol=1e-12, max_step=0.01)
    return sol.y[:, -1]
for b in (-0.9, -0.5, 0.0, 0.3, 0.8):
    y = run([-np.sqrt(1 - b * b) + 1e-12, b], [1.0, 0.0])
    assert abs(y[0]) < 1e-5 and abs(y[1] + 1) < 1e-5, (b, y)
print('OK 数値積分：平行入射はすべて (0, -a) から出る')
for ang in (10, 45, 80, 120, 170):
    u = [np.cos(np.radians(ang)), np.sin(np.radians(ang))]
    y = run([0, -1 + 1e-9], u)
    assert abs(y[2] - 1) < 1e-5 and abs(y[3]) < 1e-5, (ang, y)
print('OK 数値積分：(0, -a) からの粒子はすべて +x 向きで出る')
g = field('5')
assert g.replace(' ', '') == '(0,-a)', g
print('OK (5) 入力用の答え (0, -a)')
check_tex('6', '平行なビーム')
