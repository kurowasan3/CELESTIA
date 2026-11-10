# 検算スクリプト（problem.json の入力用の答えを sympy・numpy で独立に導いた値と照合する）
import json, re, random, math
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ESY-ROT-207'
BASE = Path(__file__).resolve().parent.parent
PDIR = BASE / 'content' / 'problems' / PID[:3] / PID
PROB = json.loads((PDIR / 'problem.json').read_text(encoding='utf-8'))
GRADING = {a['id']: a.get('grading') for a in PROB['answer']['parts']}
ANS = {a['id']: a['tex'] for a in PROB['answer']['parts']}
GREEK = {'θ': 'theta', 'ω': 'omega', 'Φ': 'Phi', 'τ': 'tau', 'Δ': 'Delta', 'φ': 'phi', 'λ': 'lam',
         'α': 'alpha', 'β': 'beta', 'σ': 'sigma', 'μ': 'mu'}
FUNCS = {'sqrt', 'sin', 'cos', 'tan', 'log', 'ln', 'exp', 'pi'}
TR = standard_transformations + (implicit_multiplication, convert_xor)


def parse(src):
    s = src.strip()
    for k, v in GREEK.items():
        s = s.replace(k, ' ' + v + ' ')
    s = s.replace('π', ' pi ')
    s = re.sub(r'\s+_', '_', s)
    s = re.sub(r'(\d)\s*([A-Za-z(])', r'\1*\2', s)
    s = re.sub(r'(?<![\w.])(\d+(?:\.\d+)?)\*10\^\(?(-?\d+)\)?', r'(\1*10**(\2))', s)
    names = set(re.findall(r'[A-Za-z][A-Za-z0-9_]*', s)) - FUNCS
    loc = {n: sp.Symbol(n, positive=True) for n in names}
    loc.update({'sqrt': sp.sqrt, 'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'log': sp.log, 'ln': sp.log,
                'exp': sp.exp, 'pi': sp.pi})
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
    names = {}
    for x in sorted(a.free_symbols | b.free_symbols, key=str):
        names.setdefault(x.name, []).append(x)
    rnd = random.Random(1)
    for _ in range(n):
        sub = {}
        for nm, xs in names.items():
            v = rnd.uniform(0.1, 1.4) if nm in ('theta', 'alpha', 'phi') else rnd.uniform(0.3, 1.7)
            for x in xs:
                sub[x] = v
        va, vb = complex(a.subs(sub).evalf()), complex(b.subs(sub).evalf())
        if abs(va - vb) > 1e-9 * max(1, abs(va), abs(vb)):
            return False
    return True


def S(names):
    return sp.symbols(names, positive=True)


def check(pid, derived, k=0):
    f = field(pid, k)
    ok = same(parse(f), derived)
    print(('OK ' if ok else 'NG ') + f'({pid}) 入力用の答え {f}  ／ 導出 {sp.simplify(derived)}')
    assert ok, f'({pid}) が一致しません'


def rnd_sig(x, sig):
    return float(f'{x:.{sig - 1}e}')


def check_num(pid, value, sig=2, k=0):
    f = field(pid, k)
    got = float(parse(f))
    want = rnd_sig(value, sig) if sig else value
    ok = abs(got - want) <= 1e-9 * max(1, abs(want))
    print(('OK ' if ok else 'NG ') + f'({pid}) 入力用の答え {f} = {got:g}  ／ 計算 {value:.6g} → {want:g}')
    assert ok, f'({pid}) の数値が一致しません'


def check_tex(pid, *needles):
    for nd in needles:
        assert nd in ANS[pid], f'({pid}) の答え {ANS[pid]} に {nd} がありません'
    print(f'OK ({pid}) 答えの文に ' + '，'.join(needles))


def fig(name):
    return json.loads((PDIR / f'{name}.fig.json').read_text(encoding='utf-8'))

# ----------------------------------------------------------------
from scipy.integrate import solve_ivp
B, w, a, R, M, g, b, v = S('B omega a R M g b v')
# 起電力：棒の各点 r の速さ rω、dV = B rω dr を積分
r = sp.Symbol('r', positive=True)
V = sp.integrate(B * r * w, (r, 0, a))
check('1', V)
I = V / R
check('2', I, 0)
N = sp.integrate(I * B * r, (r, 0, a))     # 各部分 dr の力 IB dr、腕の長さ r
check('2', N, 1)
beta = g - (N.subs(w, v / b) / b) / M
check('3', beta)
vf = sp.solve(sp.Eq(beta, 0), v)[0]
check('4', vf)
P = sp.simplify((I**2 * R).subs(w, vf / b))
check('5', P)
assert sp.simplify(P - M * g * vf) == 0
If = sp.simplify(I.subs(w, vf / b))
check('6', If)
assert not If.has(R) and sp.simplify(vf / R).has(R) is False
# 数値積分で終端速度へ近づくこと
p = {B: 0.5, a: 0.2, R: 0.1, M: 0.05, g: 9.8, b: 0.02}
fb = sp.lambdify(v, beta.subs(p))
s = solve_ivp(lambda t, y: [fb(y[0])], [0, 50], [0.0], rtol=1e-11, atol=1e-13)
assert abs(s.y[0][-1] - float(vf.subs(p))) < 1e-8
print('OK 数値：終端速度', s.y[0][-1])
print('ALL OK')
