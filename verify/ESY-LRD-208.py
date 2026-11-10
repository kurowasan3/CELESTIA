# 検算スクリプト（problem.json の入力用の答えを sympy・numpy で独立に導いた値と照合する）
import json, re, random, math
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ESY-LRD-208'
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
B, l, L, m, v0, x = S('B l L m v_0 x')
I = B * l * x / L
check('1', I)
w = sp.sqrt(B**2 * l**2 / (m * L))
check('2', w)
xm = v0 / w
check('3', xm)
Im = I.subs(x, xm)
check('4', Im)
assert sp.simplify(L * Im**2 / 2 - m * v0**2 / 2) == 0
check('5', sp.pi / (2 * w))
# 数値：回路方程式 L dI/dt = B l v と運動方程式 m dv/dt = -I B l を直接積分
p = {B: 0.8, l: 0.5, L: 0.3, m: 0.2, v0: 1.5}
Bn, ln_, Ln, mn, vn = [float(p[k]) for k in (B, l, L, m, v0)]
def f(t, y):
    xx, vv, II = y
    return [vv, -II * Bn * ln_ / mn, Bn * ln_ * vv / Ln]
ev = lambda t, y: y[1]; ev.direction = -1
s = solve_ivp(f, [0, 10], [0, vn, 0], events=ev, rtol=1e-12, atol=1e-14)
t1 = s.t_events[0][0]; xs = s.y_events[0][0]
assert abs(t1 - float((sp.pi / (2 * w)).subs(p))) < 1e-8
assert abs(xs[0] - float(xm.subs(p))) < 1e-8 and abs(xs[2] - float(Im.subs(p))) < 1e-8
assert abs(Ln * xs[2] - Bn * ln_ * xs[0]) < 1e-9
print('ALL OK')
