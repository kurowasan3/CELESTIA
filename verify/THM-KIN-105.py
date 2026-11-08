# THM-KIN-105 の答えを独立に導いて problem.json の入力用の答えと照合する（python3 -I verify/THM-KIN-105.py）
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
p, g = fields('THM-KIN-105')
n1, n2, R, T = S('n_1'), S('n_2'), S('R'), S('T')
MHe, MO = sp.Rational(4), sp.Rational(32)
check(g['1'][0], sp.sqrt(3 * R * T / MHe) / sp.sqrt(3 * R * T / MO), '(1) 二乗平均速度の比')
U = sp.Rational(3, 2) * n1 * R * T + sp.Rational(5, 2) * n2 * R * T
check(g['2'][0], U, '(2) 内部エネルギー')
CV = sp.diff(U, T) / (n1 + n2)
check(g['3'][0], CV, '(3) C_V')
# 定圧：Q = ΔU + P ΔV = ΔU + (n1+n2) R ΔT
Cp = (sp.diff(U, T) + (n1 + n2) * R) / (n1 + n2)
gam = sp.simplify(Cp / CV)
check(g['4'][0], gam, '(4) 比熱比')
x = sp.Symbol('x', positive=True)
r = sp.solve(sp.Eq(gam.subs(n2, x * n1), sp.Rational(3, 2)), x)
assert len(r) == 1
check(g['5'][0], r[0], '(5) n2/n1')
print('THM-KIN-105: すべて一致')
