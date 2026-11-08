# 検算スクリプト（problem.json の入力用の答えを sympy で独立に導いた値と照合する）
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ELM-LOR-201'
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
m, q, B, v, th = S('m q B v theta')
check('1', m * v / (q * B))
check('2', 2 * m * v * sp.cos(th) / (q * B))
check('3', (sp.pi - 2 * th) * m / (q * B))
check('5', (sp.pi + 2 * th) * m / (q * B))
# 数値シミュレーションで確かめる（B は -z）
def run(qq, mm=1.3, BB=0.7, vv=1.1, tt=0.5):
    def f(t, y):
        vx, vy = y[2], y[3]
        # F = q v × B,  B = (0,0,-BB)
        return [vx, vy, qq * (vy * (-BB)) / mm, qq * (-vx * (-BB)) / mm]
    ev = lambda t, y: y[0]; ev.terminal = True; ev.direction = -1
    s = solve_ivp(f, [0, 100], [0, 0, vv * np.cos(tt), vv * np.sin(tt)], events=ev, rtol=1e-11, atol=1e-12, max_step=0.01)
    return s.t_events[0][0], s.y_events[0][0]
sub = {m: 1.3, B: 0.7, v: 1.1, th: 0.5, q: 0.9}
t1, y1 = run(0.9)
assert abs(t1 - float(((sp.pi - 2 * th) * m / (q * B)).subs(sub))) < 1e-6
assert abs(y1[1] - float((2 * m * v * sp.cos(th) / (q * B)).subs(sub))) < 1e-6
ang = np.degrees(np.arctan2(y1[3], -y1[2]))
assert abs(ang - np.degrees(0.5)) < 1e-4 and y1[3] > 0   # x 軸の負の向きと θ、y 成分正
t2, _ = run(-0.9)
assert abs(t2 - float(((sp.pi + 2 * th) * m / (q * B)).subs(sub))) < 1e-6
print('OK 数値シミュレーション：t1, y1, 出る向き, t2 が一致')
check_tex('4', '負の向きと角 $\\theta$')
