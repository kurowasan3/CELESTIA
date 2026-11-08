# 検算スクリプト（problem.json の入力用の答えを sympy で独立に導いた値と照合する）
import json, re, random
from pathlib import Path
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication, convert_xor

PID = 'ELM-IND-202'
PROB = json.loads((Path(__file__).resolve().parent.parent / 'content' / 'problems' / PID[:3] / PID / 'problem.json').read_text(encoding='utf-8'))
GRADING = {a['id']: a.get('grading') for a in PROB['answer']['parts']}
GREEK = {'μ': 'mu', 'θ': 'theta', 'ω': 'omega', 'Φ': 'Phi', 'τ': 'tau', 'Δ': 'Delta', 'φ': 'phi'}
FUNCS = {'sqrt', 'sin', 'cos', 'tan', 'log', 'exp', 'pi'}
TR = standard_transformations + (implicit_multiplication, convert_xor)


def S(names):
    """記号（名前は入力用の答えと同じ：μ_0 → mu_0）"""
    return sp.symbols(names, positive=True)


def parse(src):
    s = src.strip()
    for k, v in GREEK.items():
        s = s.replace(k, ' ' + v + ' ' if v != 'mu' else v)
    s = s.replace('π', ' pi ')
    s = re.sub(r'\s+_', '_', s)
    s = re.sub(r'(\d)\s*([A-Za-z(])', r'\1*\2', s)
    names = set(re.findall(r'[A-Za-z][A-Za-z0-9_]*', s)) - FUNCS
    loc = {n: sp.Symbol(n, positive=True) for n in names}
    loc.update({'sqrt': sp.sqrt, 'sin': sp.sin, 'cos': sp.cos, 'tan': sp.tan, 'log': sp.log, 'exp': sp.exp, 'pi': sp.pi})
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
import numpy as np
# 双極子（N 極を下：磁気モーメント m は -z 向き）による、コイル（z=0、半径 a）を貫く磁束を数値で求め、
# 誘導起電力の向き（右ねじ：+z 向きの磁場をつくる電流 = 上から見て反時計回り）を判定する
def flux(z0, mz, a=1.0, n=400):
    # 点双極子 (0,0,z0) の磁場の z 成分をコイル面で積分
    rs = (np.arange(n) + 0.5) * a / n
    dz = -z0
    r3 = (rs**2 + dz**2) ** 1.5
    r5 = (rs**2 + dz**2) ** 2.5
    Bz = (3 * dz * (mz * dz) / r5 - mz / r3)
    return np.sum(Bz * 2 * np.pi * rs * a / n)
def emf_sign(z_from, z_to, mz):
    dPhi = flux(z_to, mz) - flux(z_from, mz)
    return -np.sign(dPhi)    # 起電力 > 0 は +z 向きの磁場をつくる向き（上から見て反時計回り）
assert flux(2.0, -1) < 0   # N 極を下：コイルを下向きに貫く
s1 = emf_sign(2.0, 1.5, -1); print('(1) N 極を近づける：', '反時計回り' if s1 > 0 else '時計回り'); assert s1 > 0
s3 = emf_sign(1.5, 2.0, +1); print('(3) S 極を下にして遠ざける：', '反時計回り' if s3 > 0 else '時計回り'); assert s3 > 0
# (4) 反時計回りの誘導電流（コイルの軸上の磁場 B_z(z) > 0、z とともに減る）と、下向きの磁気モーメント (-z) の磁石：
# エネルギー U = -m·B より F_z = d(m·B)/dz = d(-B_z)/dz
zz = sp.symbols('z', positive=True)
Bz = 1 / (1 + zz**2) ** sp.Rational(3, 2)
F_on_magnet = sp.diff(-Bz, zz).subs(zz, 1.5)
print('(4) 磁石が受ける力の z 成分（相対値）', F_on_magnet)
assert F_on_magnet > 0
check_tex('1', '反時計回り'); check_tex('2', '流れない'); check_tex('3', '反時計回り'); check_tex('4', '上向き'); check_tex('5', '大きく')
assert PROB['predict']['answer'] == 'c'
print('OK (2) 磁束が変化しない → 0、(5) 変化が速いほど大きい')
