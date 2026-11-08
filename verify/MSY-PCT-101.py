# --- 共通：problem.json の入力用の答え（grading）を sympy の式に直し，独立に導いた式と数値で照合する
import json, re, random, pathlib
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor
ROOT = pathlib.Path(__file__).resolve().parent.parent
TR = standard_transformations + (implicit_multiplication_application, convert_xor)


def grading(pid, part):
    f = next(ROOT.glob(f'content/problems/*/{pid}/problem.json'))
    p = json.loads(f.read_text(encoding='utf-8'))
    g = next(x for x in p['answer']['parts'] if x['id'] == part)['grading']
    out = []
    for fld in re.split(r'\s;\s', g):
        if '：' in fld:
            fld = fld.split('：', 1)[1]
        fld = re.sub(r'\[[^\]]*\]\s*$', '', fld).strip()
        fld = re.sub(r'倍$', '', fld).strip()
        s = fld.replace('π', 'pi').replace('θ', 'theta').replace('μ', 'mu').replace('ω', 'omega')
        s = re.sub(r'_\{([A-Za-z0-9]+)\}', r'_\1', s)
        names = set(re.findall(r'[A-Za-z_][A-Za-z_0-9]*', s)) - {'sqrt', 'sin', 'cos', 'tan', 'pi'}
        loc = {n: sp.Symbol(n, positive=True) for n in names}
        out.append(parse_expr(s, local_dict=loc, transformations=TR))
    return out


def same(a, b, trials=6):
    a, b = sp.sympify(a), sp.sympify(b)
    allsyms = a.free_symbols | b.free_symbols
    ok = 0
    for _ in range(trials * 20):
        if ok >= trials:
            break
        vals = {}
        for n in {s.name for s in allsyms}:
            if n == 'e':
                vals[n] = random.uniform(0.1, 0.9)
            elif n == 'theta':
                vals[n] = random.uniform(0.2, 1.3)
            elif n in ('mu', 'mu_0'):
                vals[n] = random.uniform(0.05, 0.3)
            else:
                vals[n] = random.uniform(1.0, 3.0)
        sub = {s: vals[s.name] for s in allsyms}
        x, y = complex(a.subs(sub).evalf()), complex(b.subs(sub).evalf())
        if abs(x.imag) > 1e-12 or abs(y.imag) > 1e-12:
            continue      # 物理的でない値の組（平方根の中が負）は使わない
        ok += 1
        if abs(x - y) > 1e-9 * max(1, abs(x), abs(y)):
            return False
    return ok >= trials


def check(pid, part, *expected):
    got = grading(pid, part)
    assert len(got) == len(expected), (pid, part, got)
    for g, e in zip(got, expected):
        assert same(g, e), f'{pid} ({part}): grading {g} ≠ 導出 {sp.simplify(e)}'
    print(f'{pid} ({part}) OK:', ', '.join(str(sp.nsimplify(e) if not e.free_symbols else sp.simplify(e)) for e in map(sp.sympify, expected)))
# ---
m, M, g, l = sp.symbols('m M g l', positive=True)
v, V = sp.symbols('v V', positive=True)
sol = sp.solve([sp.Eq(m * v, M * V), sp.Eq(m * g * l, m * v**2 / 2 + M * V**2 / 2)], [v, V], dict=True)[0]
check('MSY-PCT-101', '1', sol[v])
check('MSY-PCT-101', '2', sol[V])
S = m * g + m * (sol[v] + sol[V])**2 / l
check('MSY-PCT-101', '3', S)
check('MSY-PCT-101', '4', 2 * m * l / (M + m))
# ラグランジュ方程式で数値積分し，最下点の張力と台車の変位を確かめる
import numpy as np
mm, MM, gg, ll = 0.3, 1.0, 9.8, 0.5
# 一般化座標 X（台車）, φ（鉛直から，φ=π/2 で糸が水平に右）。小球 = (X + l sinφ, -l cosφ)
def acc(X, ph, Xd, phd):
    A = np.array([[MM + mm, mm * ll * np.cos(ph)], [mm * ll * np.cos(ph), mm * ll**2]])
    b = np.array([mm * ll * np.sin(ph) * phd**2, -mm * gg * ll * np.sin(ph)])
    return np.linalg.solve(A, b)
st = np.array([0.0, np.pi / 2, 0.0, 0.0]); dt = 1e-5; Tbot = None
def F(s):
    a = acc(*s); return np.array([s[2], s[3], a[0], a[1]])
prev = st.copy()
while True:
    k1 = F(st); k2 = F(st + dt / 2 * k1); k3 = F(st + dt / 2 * k2); k4 = F(st + dt * k3)
    new = st + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    if Tbot is None and st[1] > 0 and new[1] <= 0:
        a = acc(*new)
        # 張力：小球の運動方程式（糸の向き成分）から
        ph, phd = new[1], new[3]
        ax = a[0] + ll * (np.cos(ph) * a[1] - np.sin(ph) * phd**2); ay = ll * (np.sin(ph) * a[1] + np.cos(ph) * phd**2)
        Tbot = mm * ay + mm * gg
    if st[1] > -np.pi / 2 + 1e-3 and new[1] <= -np.pi / 2 + 1e-3:
        Xend = new[0]; break
    st = new
assert abs(Tbot - float(S.subs({m: mm, M: MM, g: gg, l: ll}))) < 1e-4, Tbot
assert abs(abs(Xend) - 2 * mm * ll / (MM + mm)) < 1e-3, Xend
