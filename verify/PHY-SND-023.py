# PHY-SND-023 : ピストンが送り出す圧縮波から音速を導き，混合気体の組成を定常波の節の間隔から決める
import sympy as sp

g, p0, rho, wv, c, t, S = sp.symbols('gamma p_0 rho w c t S', positive=True)
# (1) 断熱：p V^g 一定，体積減少率 w/c
h = sp.symbols('h')
dp_adiab = sp.series(p0*(1 - h)**(-g) - p0, h, 0, 2).removeO().subs(h, wv/c)
print('dp adiabatic:', sp.simplify(dp_adiab))
# (2) 運動量：Δp S t = (rho S c t) w
dp_mom = rho*c*wv
csol = sp.solve(sp.Eq(dp_adiab, dp_mom), c)
print('c =', csol)
assert sp.simplify(csol[0]**2 - g*p0/rho) == 0

# (5) 混合気体  He: x, N2: 1-x
x, R, T = sp.symbols('x R T', positive=True)
Cv = sp.Rational(3, 2)*R*x + sp.Rational(5, 2)*R*(1-x)
gam = (Cv + R)/Cv
M = 4*x + 28*(1-x)
c2 = gam*R*T/M
c2N = sp.Rational(7, 5)*R*T/28
sols = sp.solve(sp.Eq(c2/c2N, sp.Rational(15, 8)), x)
print('x solutions:', sols)
assert sp.Rational(1, 2) in sols
# 単調性（0<=x<=1 で比が単調増加 → 解は一意）
ratio = sp.simplify(c2/c2N)
dr = sp.diff(ratio, x)
import numpy as np
xs = np.linspace(0, 1, 1001)
vals = [float(dr.subs(x, xv)) for xv in xs]
assert min(vals) > 0
print('ratio at x=1:', float(ratio.subs(x, 1)))
# (6) 等温と誤解した場合
xi = sp.solve(sp.Eq(28/M, sp.Rational(15, 8)), x)
print('isothermal inference x =', xi, [float(v) for v in xi])
assert xi == [sp.Rational(49, 90)]
# 数値：f=500 Hz, T=300 K の N2 の節の間隔
cN = float(sp.sqrt(sp.Rational(7, 5)*8.31*300/0.028))
print('c_N2 =', cN, ' d0 =', cN/2/500)
# (4) 壁の位置は変位の節＝圧力の腹．最も近い圧力の節は λ/4 = d/2
# 定常波  ξ = A sin(kx') cos(ωt) (x' は壁からの距離, 壁で ξ=0)，p' = -K ∂ξ/∂x ∝ cos(kx')
k = sp.symbols('k', positive=True); xp = sp.symbols('xp', nonnegative=True)
pp = sp.diff(sp.sin(k*xp), xp)
print('pressure amplitude ∝', pp, ' first zero at kx=', sp.solve(sp.cos(k*xp), xp))
print('ALL OK')
