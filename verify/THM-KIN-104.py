# THM-KIN-104 の答えを独立に導いて problem.json の入力用の答えと照合する（python3 -I verify/THM-KIN-104.py）
# --- 共通：problem.json の入力用の答え（grading）を sympy の式に読み直して照合する
import json, re, sys
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
GREEK = {'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'ρ': 'rho', 'Δ': 'Delta', 'θ': 'theta', 'π': 'pi'}
TOK = re.compile(r'\s*(sqrt|\d+(?:\.\d+)?|[A-Za-zα-ωΔ](?:_(?:\{[^}]*\}|[A-Za-z0-9]+))?|\*\*|[-+*/^()])')

def to_sympy(src):
    s = src.strip()
    s = re.sub(r'\s*\[[^\]]*\]\s*$', '', s)      # 単位
    s = re.sub(r'倍$', '', s).strip()
    toks, pos = [], 0
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m: raise ValueError(f'読めない：{s[pos:]}')
        toks.append(m.group(1)); pos = m.end()
    out, names = [], {}
    def is_val(t): return t not in ('sqrt', '**', '^') and (t == ')' or re.match(r'[\dA-Za-zα-ωΔ]', t))
    def starts_val(t): return t in ('(', 'sqrt') or re.match(r'[\dA-Za-zα-ωΔ]', t)
    prev = None
    for t in toks:
        if prev is not None and is_val(prev) and starts_val(t): out.append('*')
        if t == '^': out.append('**')
        elif re.match(r'[A-Za-zα-ωΔ]', t) and t != 'sqrt':
            base = t.replace('{', '').replace('}', '')
            py = ''.join(GREEK.get(c, c) for c in base)
            names[py] = sp.Symbol(base.replace('_', '_'), positive=True)
            if py == 'pi': names[py] = sp.pi
            out.append(py)
        else: out.append(t)
        prev = t
    return sp.sympify(''.join(out), locals={**names, 'sqrt': sp.sqrt}, rational=True)

def fields(pid):
    p = json.loads((ROOT / 'content/problems/THM' / pid / 'problem.json').read_text(encoding='utf-8'))
    res = {}
    for a in p['answer']['parts']:
        g = a.get('grading')
        if not g: continue
        res[a['id']] = [to_sympy(re.sub(r'^[^：]{1,80}：', '', x.strip())) for x in g.split(' ; ')]
    return p, res

def S(name): return sp.Symbol(name, positive=True)

def check(got, want, label):
    got, want = sp.sympify(got), sp.sympify(want)
    d = sp.simplify(got - want)
    ok = d == 0
    if not ok and d.free_symbols:   # 平方根の整理の違いなど：正の乱数を何組か代入して確かめる
        import random
        rnd = random.Random(0)
        ok = all(abs(complex(d.subs({s_: rnd.uniform(0.5, 3) for s_ in d.free_symbols}).evalf())) < 1e-9 for _ in range(20))
    print(('OK  ' if ok else 'NG  ') + f'{label}: problem.json = {got}   derived = {sp.simplify(want)}')
    if not ok: sys.exit(f'食い違い：{label}')

# --- この問題
p, g = fields('THM-KIN-104')
import numpy as np
m, L, vx, N, w, A, k, T, R = [S(x) for x in 'm L v_x N w A k T R'.split()]
# (1) 力積 2 m vx × 回数 vx/(2L)
check(g['1'][0], 2 * m * vx * vx / (2 * L), '(1) 1個の分子の平均の力')
# (2) 平面内で等方的：<vx^2> = <v^2>/2 を数値で確かめてから
rng = np.random.default_rng(1)
th = rng.uniform(0, 2 * np.pi, 400000); sp_ = rng.rayleigh(1.0, 400000)
ratio = np.mean((sp_ * np.cos(th))**2) / np.mean(sp_**2)
assert abs(ratio - 0.5) < 0.005, ratio
f = N * m * (w / 2) / L / L
check(g['2'][0], f.subs(L, sp.sqrt(A)), '(2) f')
# (3) fA = NkT → (1/2) m w = ?
Ek = sp.solve(sp.Eq(f.subs(L, sp.sqrt(A)) * A, N * k * T), w)[0] * m / 2
check(g['3'][0], Ek, '(3) 平均運動エネルギー')
# (4)(5) 1 mol：U = N_A * Ek = R T
dT = sp.Symbol('dT')
U = lambda TT: R * TT
CA = sp.simplify((U(T + dT) - U(T)) / dT)
check(g['4'][0], CA, '(4) C_A')
# f 一定：W = f ΔA = R ΔT
Cf = sp.simplify((U(T + dT) - U(T) + R * dT) / dT)
check(g['5'][0], Cf, '(5) C_f')
check(g['5'][1], Cf / CA, '(5) 比')
print('THM-KIN-104: すべて一致')
