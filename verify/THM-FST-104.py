# THM-FST-104 の答えを独立に導いて problem.json の入力用の答えと照合する（python3 -I verify/THM-FST-104.py）
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
p, g = fields('THM-FST-104')
P0, V0, T0 = S('P_0'), S('V_0'), S('T_0')
n_R = P0 * V0 / T0                       # nR
V = sp.Symbol('V', positive=True)
Pline = sp.expand(3 * P0 + (P0 - 3 * P0) / (2 * V0) * (V - V0))     # B(V0, 3P0) と C(3V0, P0) を通る直線
# (1) 正味の仕事：各過程の ∫P dV の和
W_BC = sp.integrate(Pline, (V, V0, 3 * V0))
W_CA = sp.integrate(P0, (V, 3 * V0, V0))
W = sp.simplify(W_BC + W_CA)
check(g['1'][0], W, '(1) W')
# (2) 温度最大
PV = Pline * V
Vmax = sp.solve(sp.diff(PV, V), V)[0]
check(g['2'][0], Vmax, '(2) 最高温度の体積')
check(g['2'][1], sp.simplify(PV.subs(V, Vmax) / n_R), '(2) 最高温度')
# (3) dQ/dV = (3/2) d(PV)/dV + P
dQdV = sp.Rational(3, 2) * sp.diff(PV, V) + Pline
Vsw = sp.solve(dQdV, V)[0]
check(g['3'][0], Vsw, '(3) 吸熱→放熱の切り替わり')
# (4) 吸熱の和：dQ>0 の部分だけを積分
Q_AB = sp.Rational(3, 2) * (3 * P0 * V0 - P0 * V0)
Q_BC_in = sp.integrate(dQdV, (V, V0, Vsw))
Q_BC_out = sp.integrate(dQdV, (V, Vsw, 3 * V0))
Q_CA = sp.Rational(5, 2) * P0 * (V0 - 3 * V0)
assert Q_AB > 0 and Q_CA < 0 and sp.simplify(Q_BC_out / (P0 * V0)) < 0
Qin = sp.simplify(Q_AB + Q_BC_in)
check(g['4'][0], Qin, '(4) Q_in')
assert sp.simplify(Q_AB + Q_BC_in + Q_BC_out + Q_CA - W) == 0      # 1サイクルで ΔU = 0
check(g['5'][0], sp.simplify(W / Qin), '(5) 熱効率')
print('THM-FST-104: すべて一致')
