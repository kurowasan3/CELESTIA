# THM-KIN-205 の答えを独立に導いて problem.json の入力用の答えと照合する（python3 -I verify/THM-KIN-205.py）
# --- 共通：problem.json の入力用の答え（grading）を sympy の式に読み直して照合する
import json, re, sys
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
GREEK = {'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'ρ': 'rho', 'Δ': 'Delta', 'θ': 'theta', 'π': 'pi', 'μ': 'mu'}
FUNCS = {'sqrt': sp.sqrt, 'log': sp.log, 'cos': sp.cos, 'sin': sp.sin}
TOK = re.compile(r'\s*(sqrt|log|cos|sin|\d+(?:\.\d+)?|[A-Za-zα-ωΔ](?:_(?:\{[^}]*\}|[A-Za-z0-9]+))?|\*\*|[-+*/^()])')


def to_sympy(src):
    s = src.strip()
    s = re.sub(r'\s*\[[^\]]*\]\s*$', '', s)      # 単位
    s = re.sub(r'倍$', '', s).strip()
    toks, pos = [], 0
    while pos < len(s):
        m = TOK.match(s, pos)
        if not m:
            raise ValueError(f'読めない：{s[pos:]}')
        toks.append(m.group(1)); pos = m.end()
    out, names = [], {}
    def is_val(t): return t not in FUNCS and t not in ('**', '^') and (t == ')' or re.match(r'[\dA-Za-zα-ωΔ]', t))
    def starts_val(t): return t == '(' or t in FUNCS or re.match(r'[\dA-Za-zα-ωΔ]', t)
    prev = None
    for t in toks:
        if prev is not None and is_val(prev) and starts_val(t):
            out.append('*')
        if t == '^':
            out.append('**')
        elif t in FUNCS:
            out.append(t)
        elif re.match(r'[A-Za-zα-ωΔ]', t):
            base = t.replace('{', '').replace('}', '')
            py = ''.join(GREEK.get(c, c) for c in base)
            names[py] = sp.pi if py == 'pi' else sp.Symbol(base, positive=True)
            out.append(py)
        else:
            out.append(t)
        prev = t
    return sp.sympify(''.join(out), locals={**names, **FUNCS}, rational=True)


def fields(pid):
    p = json.loads((ROOT / 'content/problems/THM' / pid / 'problem.json').read_text(encoding='utf-8'))
    res = {}
    for a in p['answer']['parts']:
        g = a.get('grading')
        if not g:
            continue
        res[a['id']] = [to_sympy(re.sub(r'^[^：]{1,80}：', '', x.strip())) for x in g.split(' ; ')]
    return p, res


def S(name):
    return sp.Symbol(name, positive=True)


def check(got, want, label, rel=None):
    """rel を与えると相対誤差 rel 以内（有効数字で丸めた値の照合）"""
    got, want = sp.sympify(got), sp.sympify(want)
    if rel is not None:
        g, w = float(got), float(want)
        ok = abs(g - w) <= rel * abs(w)
    else:
        d = sp.simplify(got - want)
        ok = d == 0
        if not ok and d.free_symbols:
            import random
            rnd = random.Random(0)
            ok = all(abs(complex(d.subs({s_: rnd.uniform(0.5, 3) for s_ in d.free_symbols}).evalf())) < 1e-9 for _ in range(20))
        elif not ok:
            ok = abs(complex(d.evalf())) < 1e-12
    print(('OK  ' if ok else 'NG  ') + f'{label}: problem.json = {got}   derived = {sp.N(want, 8) if rel is not None else sp.simplify(want)}')
    if not ok:
        sys.exit(f'食い違い：{label}')


def sig(x, n=2):
    """有効数字 n 桁に丸めた値（問題文で有効数字を指定した答えの照合用）"""
    x = float(x)
    from math import floor, log10
    k = floor(log10(abs(x)))
    return round(x, -k + n - 1)


def check_sig(got, want, label, n=2):
    """problem.json の値が，導いた値を有効数字 n 桁に丸めたものと一致するか"""
    r = sig(want, n)
    ok = abs(float(got) - r) <= 1e-9 * abs(r)
    print(('OK  ' if ok else 'NG  ') + f'{label}: problem.json = {got}   derived = {float(want):.6g} → 有効数字{n}桁 {r:g}')
    if not ok:
        sys.exit(f'食い違い：{label}')

# --- この問題
p, g = fields('THM-KIN-205')
import numpy as np
n, m, v, S_, th = S('n'), S('m'), S('v'), S('S'), sp.Symbol('θ', positive=True)
check(g['1'][0], n * v * S_, '(1) 単位時間に当たる数')
check(g['2'][0], (2 * m * v) * n * v * S_, '(2) 弾性衝突の力')
check(g['3'][0], (m * v) * n * v * S_, '(3) 付着の力')
# (4) ベクトルで独立に計算：分子線 +x，板の法線（来る側）は -x から θ 回転
def F4(theta):
    nhat = np.array([-np.cos(theta), np.sin(theta)])       # 来る側を向く単位法線
    vin = np.array([1.0, 0.0])
    vout = vin - 2 * vin.dot(nhat) * nhat                    # 鏡面反射
    rate = abs(vin.dot(nhat))                                # 単位面積の板が分子線をさえぎる割合（×n v S）
    dp = vout - vin                                          # 分子の運動量変化（÷m）
    return -dp * rate                                        # 板が受ける力（÷ n m v^2 S）
for tv in (0.3, 0.7, 1.1):
    f = F4(tv)
    val = float(g['4'][0].subs({n: 1, m: 1, v: 1, S_: 1, th: tv}))
    assert abs(np.linalg.norm(f) - val) < 1e-12, (tv, f, val)
    nhat = np.array([-np.cos(tv), np.sin(tv)])
    assert abs(abs(f.dot(nhat)) - np.linalg.norm(f)) < 1e-12 and f.dot(nhat) < 0   # 法線の逆向き
print('OK  (4) 斜めの板：大きさ 2nmv^2S cos^2θ，向きは法線の逆（数値で確認）')
# 図の向き：θ=30° で反射 120°，力 330°
f = F4(np.radians(30)); vout = np.array([1, 0]) - 2 * np.array([1, 0]).dot(np.array([-np.cos(np.radians(30)), np.sin(np.radians(30))])) * np.array([-np.cos(np.radians(30)), np.sin(np.radians(30))])
assert abs(np.degrees(np.arctan2(vout[1], vout[0])) - 120) < 1e-9 and abs(np.degrees(np.arctan2(f[1], f[0])) % 360 - 330) < 1e-9
# (5) 等方的な気体：P = n m v^2/3
check(g['5'][0], (n * m * v**2 / 3 * S_) / (2 * n * m * v**2 * S_), '(5) 気体の力 / F1')
print('THM-KIN-205: すべて一致')
