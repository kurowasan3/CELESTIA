# ELM-FLD-304: 答えを独立に導いて problem.json の grading と照合する（python3 -I verify/ELM-FLD-304.py）
# --- 共通：problem.json の grading を読み，sympy の式にして独立に導いた答えと照合する
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor
PID = Path(__file__).stem
P = json.loads((Path(__file__).resolve().parent.parent / 'content/problems/ELM' / PID / 'problem.json').read_text(encoding='utf-8'))
GR = {x['id']: x.get('grading') for x in P['answer']['parts']}
TR = standard_transformations + (implicit_multiplication_application, convert_xor)
GREEK = {'π': 'pi', 'ε_0': 'eps0', 'ε_r': 'epsr', 'θ_0': 'theta0', 'θ': 'theta'}
def to_expr(s):
    s = re.sub(r'\s*\[[^\[\]]+\]\s*$', '', s.strip())
    s = re.sub(r'倍$|回$', '', s).strip()
    for g, n in GREEK.items():
        s = s.replace(g, n)
    s = re.sub(r'_\{(\w+)\}', r'_\1', s)
    names = set(re.findall(r'[A-Za-z][A-Za-z0-9_]*', s)) - {'sqrt', 'pi'}
    loc = {n: sp.Symbol(n, positive=True) for n in names}
    return parse_expr(s, local_dict=loc, transformations=TR)
def fields(pid):
    out = {}
    for i, f in enumerate([x.strip() for x in GR[pid].split(' ; ')]):
        lab, _, body = f.rpartition('：')
        out[lab or i] = to_expr(body)
    return out
def same(a, b, rel=1e-9):
    a, b = sp.sympify(a), sp.sympify(b)
    syms = sorted(a.free_symbols | b.free_symbols, key=str)
    for _ in range(5):
        sub = {x: sp.Rational(random.randint(11, 97), random.randint(5, 13)) for x in syms}
        va, vb = sp.N(a.subs(sub), 30), sp.N(b.subs(sub), 30)
        if abs(va - vb) > rel * max(1, abs(vb)):
            return False
    return True
def check(pid, expected, label=None, sig=None):
    """sig を与えると，導いた値を有効数字 sig 桁に丸めたものと grading が一致するか確かめる"""
    fs = fields(pid)
    got = fs[label] if label is not None else list(fs.values())[0]
    if sig is not None:
        x = float(sp.N(expected)); expected = sp.Float(f'{x:.{sig - 1}e}')
    assert same(got, expected), f'({pid}{":" + str(label) if label is not None else ""}) grading {got} != derived {sp.N(expected)}'
    print(f'  ({pid}{":" + str(label) if label is not None else ""}) {got}  ✓')
S = lambda n: sp.Symbol(n, positive=True)
R_ = sp.Rational
# ---
# 自動採点の小問はない（すべて記述）．簡単な電荷の勘定で (1)〜(4) の符号を確かめる
assert all(g is None for g in GR.values())
ind = 1.0                     # 棒が金属板に引き寄せる電子の量（正の棒の強さに比例）
plate, leaf = -ind, +ind      # (1)
assert plate < 0 and leaf > 0
leaf = 0.0; total = plate + leaf   # (2) 接地：はくの正電荷が大地からの電子で中和
assert total < 0
leaf = total / 2              # (3) 棒を遠ざけると全体（負）に広がる（はくに半分程度）
assert leaf < 0
# (4) 正の棒を近づける強さ s：はくの電荷 = total/2 + s。s を増やすと 0 を通って正へ
leaves = [total / 2 + s for s in (0, 0.25, 0.5, 0.75, 1.0)]
assert leaves[0] < 0 and abs(leaves[2]) < 1e-12 and leaves[-1] > 0
print('OK', PID)
