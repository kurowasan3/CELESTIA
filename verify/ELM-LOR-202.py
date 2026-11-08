# 検算スクリプト（problem.json の入力用の答えを sympy で独立に導いた値と照合する）
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ELM-LOR-202'
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
m, q, B, E = S('m q B E')
u = sp.solve(sp.Eq(q * E - q * sp.Symbol('u') * B, 0), sp.Symbol('u'))[0]   # F_y = qE + q(u x̂ × B ẑ)_y = qE - quB
check('1', u)
w = q * B / m
R = u / w
check('2', R)
check('3', 2 * R)
check('4', 2 * u)
assert sp.simplify(q * E * 2 * R - m * (2 * u) ** 2 / 2) == 0
check('5', u * 2 * sp.pi / w)
check('6', u)
# 数値：運動方程式を直接解く
mm, qq, BB, EE = 1.2, 0.8, 0.9, 0.5
def f(t, y):
    vx, vy = y[2], y[3]
    return [vx, vy, qq * (vy * BB) / mm, qq * (EE - vx * BB) / mm]
T = 2 * np.pi * mm / (qq * BB)
s = solve_ivp(f, [0, T], [0, 0, 0, 0], rtol=1e-11, atol=1e-13, max_step=T / 2000, dense_output=True)
ts = np.linspace(0, T, 4001); Y = s.sol(ts)
sub = {m: mm, q: qq, B: BB, E: EE}
assert abs(Y[1].max() - float((2 * R).subs(sub))) < 1e-6
assert abs(np.hypot(Y[2], Y[3]).max() - float((2 * u).subs(sub))) < 1e-6
assert abs(Y[0][-1] - float((u * 2 * sp.pi / w).subs(sub))) < 1e-6 and abs(Y[1][-1]) < 1e-8
print('OK 数値シミュレーション：y_max, v_max, x_1 が一致')
