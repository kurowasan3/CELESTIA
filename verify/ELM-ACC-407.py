# 検算スクリプト（problem.json の入力用の答えを sympy で独立に導いた値と照合する）
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ELM-ACC-407'
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
V0, R, w, t, C, T = S('V_0 R omega t C T')
Tp = 2 * sp.pi / w
vR = V0 * sp.sin(w * t)
P = sp.integrate(vR**2 / R, (t, 0, sp.pi / w)) / Tp
check('1', sp.simplify(P))
Ve = sp.sqrt(P * R)
check('2', sp.simplify(Ve))
Iavg = sp.integrate(vR / R, (t, 0, sp.pi / w)) / Tp
check('3', sp.simplify(Iavg))
check('4', V0 * T / (R * C))
Cmin = sp.solve(sp.Eq(10 * sp.Rational(1, 50) / (1000 * C), sp.Rational(1, 2)), C)[0]
check('5', Cmin * 10**6)
dV, x = sp.symbols('Delta V', positive=True), sp.symbols('x', positive=True)
dV = sp.Symbol('Delta', positive=True) * sp.Symbol('V', positive=True)
# V0 cos x = V0 - ΔV の小さい解 ≈ sqrt(2ΔV/V0)
frac = sp.sqrt(2 * dV / V0) / (2 * sp.pi)
check('6', frac)
# 数値：ΔV/V0 = 0.05 のとき、厳密な acos と比べる
import math
exact = math.acos(1 - 0.05) / (2 * math.pi)
approx = math.sqrt(0.1) / (2 * math.pi)
assert abs(approx - 0.050) < 0.001 and abs(exact - approx) / exact < 0.02, (exact, approx)
print('OK 数値：導通時間の割合', approx, '（厳密', exact, '）')
# 数値シミュレーションで平滑回路のリプルを確かめる（V0=10, f=50, R=1k, C=400μ）
import numpy as np
V0n, fn, Rn, Cn = 10.0, 50.0, 1000.0, 400e-6
dt = 1e-6; vc = 0.0; ts = np.arange(0, 0.2, dt); hist = []
for tt in ts:
    vs = V0n * math.sin(2 * math.pi * fn * tt)
    vc -= vc / (Rn * Cn) * dt
    if vs > vc:
        vc = vs
    hist.append(vc)
h = np.array(hist[-int(0.02 / dt):])
ripple = h.max() - h.min()
assert 0.40 < ripple <= 0.50, ripple
print('OK 数値：リプル', ripple, 'V（≦ 0.50 V）')
