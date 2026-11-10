# 検算スクリプト（problem.json の入力用の答えを sympy で独立に導いた値と照合する）
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ELM-ACC-408'
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
V, R, L, C, w = S('V R L C omega')
I = V / sp.sqrt(R**2 + (w * L - 1 / (w * C))**2)
check('1', I)
w0 = 1 / sp.sqrt(L * C)
check('2', w0, k=0)
check('2', (I**2 * R).subs(w, w0), k=1)
sols = sp.solve(sp.Eq((w * L - 1 / (w * C))**2, R**2), w)
pos = [s for s in sols if sp.simplify(s.subs({R: 1, L: 1, C: 1})) > 0]
assert len(pos) == 2
w1, w2 = sorted(pos, key=lambda s: float(s.subs({R: 1, L: 1, C: 1})))
check('3', sp.simplify(w2 - w1))
assert sp.simplify(w1 * w2 - w0**2) == 0
ratio = (V / R) / (w0 * C) / V
check('4', ratio)
assert sp.simplify(ratio - w0 / (w2 - w1)) == 0
VC = (ratio * V).subs({L: sp.Rational(1, 100), C: sp.Rational(1, 10**7), R: 10, V: 1})
assert abs(float(VC) - 31.62) < 0.01
check('5', 32)
t = sp.symbols('t')
Irms = V / R
i = sp.sqrt(2) * Irms * sp.sin(w0 * t)
q = -sp.sqrt(2) * Irms / w0 * sp.cos(w0 * t)
assert sp.simplify(sp.diff(q, t) - i) == 0
U = sp.simplify(L * i**2 / 2 + q**2 / (2 * C))
assert sp.diff(U, t).simplify() == 0
check('6', U, k=0)
W = Irms**2 * R * 2 * sp.pi / w0
check('6', sp.simplify(2 * sp.pi * U / W), k=1)
