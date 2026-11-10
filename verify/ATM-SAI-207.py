# 検算スクリプト（problem.json の入力用の答えを sympy・numpy で独立に導いた値と照合する）
import json, re, random, math
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ATM-SAI-207'
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
I, S_, lam, h, c, th = S('I S lam h c theta')
n = I * S_ / (h * c / lam)
check('1', n)
F0 = n * 2 * h / lam
check('2', F0)
# 傾けた帆：光子の運動量ベクトルを鏡面反射させて変化を計算
th_ = sp.Symbol('theta', positive=True)
pin = sp.Matrix([1, 0]) * h / lam
nrm = sp.Matrix([sp.cos(th_), sp.sin(th_)])
pout = pin - 2 * pin.dot(nrm) * nrm
dp = pin - pout                               # 帆が受け取る運動量
Fv = n * sp.cos(th_) * dp
Fmag = sp.simplify(sp.sqrt(Fv.dot(Fv)))
check('3', 2 * I * S_ * sp.cos(th)**2 / c)
assert same(Fmag.subs(th_, th), 2 * I * S_ * sp.cos(th)**2 / c)
assert sp.simplify(Fv[0] * nrm[1] - Fv[1] * nrm[0]) == 0      # 力は法線方向
In, Sn, cn, g0 = 1.4e3, 1.0e4, 3.0e8, 5.9e-3
F0n = 2 * In * Sn / cn
check_num('4', F0n, 2, 0); check_num('4', F0n / 50, 2, 1)
check_num('5', 2 * In / (cn * g0), 2)
# (6) 光線に垂直な成分 cos^2 sin の最大
s = sp.Symbol('s', positive=True)
sol = [x for x in sp.solve(sp.diff(s - s**3, s), s) if x.is_positive][0]
check('6', sol)
import numpy as np
tt = np.linspace(0, np.pi / 2, 200001)
assert abs(np.sin(tt[np.argmax(np.cos(tt)**2 * np.sin(tt))]) - 1 / math.sqrt(3)) < 1e-4
print('ALL OK')
