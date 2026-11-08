# THM-KIN-102 の答えを独立に導いて problem.json の入力用の答えと照合する（python3 -I verify/THM-KIN-102.py）
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
p, g = fields('THM-KIN-102')
N, m, V, w, P = sp.symbols('N m V w P', positive=True)
rho = S('ρ')
# P = N m w/(3V)、ρ = N m / V から w を P と ρ で
w_expr = sp.solve(sp.Eq(P, N * m * w / (3 * V)), w)[0].subs(N * m, rho * V)
w_expr = sp.simplify(w_expr.subs(N, rho * V / m))
check(g['1'][0].subs(S('P'), P), w_expr, '(1) v^2 の平均')
v2 = sp.sqrt(w_expr.subs({P: 10**5, rho: sp.Rational(12, 10)}))
check(g['2'][0], v2, '(2) 二乗平均速度 [m/s]')
# 体積一定（ρ 一定）で T を 4 倍 → P も 4 倍
v3 = sp.sqrt(w_expr.subs({P: 4 * 10**5, rho: sp.Rational(12, 10)}))
check(g['3'][0], v3, '(3) 加熱後の二乗平均速度 [m/s]')
print('THM-KIN-102: すべて一致')
