# 検算スクリプト（problem.json の入力用の答えを sympy で独立に導いた値と照合する）
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ELM-ACC-204'
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
L, C, V0 = S('L C V_0')
t = sp.symbols('t', positive=True)
q = sp.Function('q')
sol = sp.dsolve(L * q(t).diff(t, 2) + q(t) / C, q(t), ics={q(0): C * V0, q(t).diff(t).subs(t, 0): 0}).rhs
w = 1 / sp.sqrt(L * C)
assert sp.simplify(sol - C * V0 * sp.cos(w * t)) == 0
check('1', 2 * sp.pi / w)
i = -sp.diff(sol, t)
check('2', sp.simplify(i.subs(t, sp.pi / (2 * w))))
t1 = sp.solve(sp.Eq(sp.cos(w * t), sp.Rational(1, 2)), t)
t1 = min(t1, key=lambda e: float(e.subs({L: 1, C: 1})))
check('3', t1)
check('4', sp.simplify(i.subs(t, t1)))
# (5)(6) 電流最大の瞬間に C2 を並列：合成 2C で振動。初期条件 q=0、i=I_max
Imax = V0 * sp.sqrt(C / L)
check('5', 2 * sp.pi * sp.sqrt(L * 2 * C))
Vp = sp.solve(sp.Eq(sp.Rational(1, 2) * (2 * C) * sp.Symbol('Vp')**2, sp.Rational(1, 2) * L * Imax**2), sp.Symbol('Vp'))
Vp = [x for x in Vp if x.subs({V0: 1, L: 1, C: 1}) > 0][0]
check('6', Vp)
# 数値：途中でつなぐ全体を積分
LL, CC, VV = 0.3, 0.02, 2.0
tA = np.pi / 2 * np.sqrt(LL * CC)
s1 = solve_ivp(lambda t, y: [y[1], -y[0] / (LL * CC)], [0, tA], [CC * VV, 0], rtol=1e-12, atol=1e-14)
qA, dqA = s1.y[0][-1], s1.y[1][-1]
s2 = solve_ivp(lambda t, y: [y[1], -y[0] / (LL * 2 * CC)], [0, 10], [qA, dqA], rtol=1e-12, atol=1e-14, max_step=1e-3)
Vmax = np.abs(s2.y[0]).max() / (2 * CC)
assert abs(Vmax - VV / np.sqrt(2)) < 1e-5 and abs(qA) < 1e-9
print('OK 数値：つないだ後の電圧の最大値', Vmax)
