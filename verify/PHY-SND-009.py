# K20103-1 検証：気柱の共鳴と加熱の仕方（定圧・定積・断熱圧縮）
import sympy as sp
import numpy as np

lA, f, TA, VA, T, ell, alpha = sp.symbols('ell_A f T_A V_A T ell alpha', positive=True)

# 振動板 x=0 が変位の腹, ふた x=ell が変位の節: ell = (2n-1) lambda/4
# 状態A: 腹が5か所 → 2n-1=9
lamA = 4*lA/9
VA_ = f*lamA
print('lambda_A =', lamA, ' V_A =', VA_)
# 密度の変動が最大の位置 = 変位の節: x = (2j-1) lambda/4 = (2j-1) ell_A/9
print('density-max positions:', [sp.Rational(2*j-1, 9) for j in range(1, 6)], '× ell_A')

# 音速 V = V_A sqrt(T/T_A)。 次数 N=2n-1 = 4 ell f / V
def N(ell_, T_):
    return 4*ell_*f/(VA_*sp.sqrt(T_/TA))
# (2) 定圧: ell = ell_A T/T_A
TB = sp.solve(sp.Eq(N(lA*T/TA, T), 11), T)[0]
print('isobaric next: T_B/T_A =', sp.simplify(TB/TA), ' ell_B/ell_A =', sp.simplify(TB/TA))
# (3) 定積: ell=ell_A
TC = sp.solve(sp.Eq(N(lA, T), 7), T)[0]
print('isochoric next: T_C/T_A =', sp.simplify(TC/TA))
# (4) 断熱圧縮 gamma=7/5: T ell^(gamma-1) = const
g = sp.Rational(7, 5)
x = sp.symbols('x', positive=True)  # x = ell/ell_A
Tad = TA*x**(-(g-1))
expr = sp.simplify(N(lA*x, Tad))
print('adiabatic N(x) =', expr)
xD = sp.solve(sp.Eq(expr, 7), x)[0]
print('adiabatic: ell_D/ell_A =', xD, '=', float(xD), ' T_D/T_A =', sp.nsimplify(Tad.subs(x, xD)/TA), float(Tad.subs(x, xD)/TA))
# (5) 新設問：N は ell と p だけで決まる。N = 9 sqrt((ell/ell_A)(p_A/p))
p_, pA = sp.symbols('p p_A', positive=True)
T_of = TA*p_*ell/(pA*lA)                 # 状態方程式 p*ell/T 一定
N_lp = sp.simplify(N(ell, T_of))
print('N(ell,p) =', N_lp)
assert sp.simplify(N_lp - 9*sp.sqrt(ell*pA/(lA*p_))) == 0
assert not N_lp.has(T) and not N_lp.has(TA)
# (2)(3)(4) をこの式で再現
assert sp.simplify(N_lp.subs({p_: pA, ell: lA*sp.Rational(121, 81)}) - 11) == 0          # 定圧
assert sp.simplify(N_lp.subs({ell: lA, p_: pA*sp.Rational(81, 49)}) - 7) == 0            # 定積: p ∝ T
xD_ = (sp.Rational(7, 9))**sp.Rational(5, 6)
assert abs(float(N_lp.subs({ell: lA*xD_, p_: pA*xD_**(-sp.Rational(7, 5))})) - 7) < 1e-12  # 断熱: p ∝ ell^-gamma
# 共鳴が続く条件 p ∝ ell, T ∝ ell^2
assert sp.simplify(N_lp.subs(p_, pA*ell/lA) - 9) == 0
print('T along p∝ell:', sp.simplify(T_of.subs(p_, pA*ell/lA)))
assert sp.simplify(T_of.subs(p_, pA*ell/lA) - TA*(ell/lA)**2) == 0
# 数値：p∝ell の経路で ell を 1→2 倍に変えても N は 9 のまま
for xv in np.linspace(1, 2, 11):
    assert abs(float(N_lp.subs({ell: lA*xv, p_: pA*xv})) - 9) < 1e-12

# 数値: 温度を少しずつ上げて共鳴（N が奇数整数）となる温度を走査
fn, lAn, TAn = 1000.0, 1.0, 300.0
VAn = fn*4*lAn/9
Ts = np.linspace(300, 500, 2000001)
Nn = 4*(lAn*Ts/TAn)*fn/(VAn*np.sqrt(Ts/TAn))
idx = np.where(np.diff(np.floor((Nn-1)/2)) > 0)[0]
print('isobaric scan first resonance after A at T =', Ts[idx[0]+1], ' formula', TAn*121/81)
Ts2 = np.linspace(300.01, 600, 3000001)
Nn2 = 4*lAn*fn/(VAn*np.sqrt(Ts2/TAn))
idx2 = np.where(np.diff(np.floor((Nn2-1)/2)) < 0)[0]
print('isochoric scan first resonance after A at T =', Ts2[idx2[0]], ' formula', TAn*81/49)
