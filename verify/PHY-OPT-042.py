# PHY-OPT-042 検証：点光源と目が軸上にあるときのCDの虹色の輪
import sympy as sp, numpy as np
from scipy.optimize import brentq
# (1) 基板内：n d (sin an + sin bn) と sin a = n sin an
n,d,an,bn=sp.symbols('n d alpha_n beta_n',positive=True)
print('(1)',sp.simplify(n*d*(sp.sin(an)+sp.sin(bn)) - d*(n*sp.sin(an)+n*sp.sin(bn))))
D=1.60e-3 # mm
def S(r,h): return r/np.hypot(r,h)
def ring(m,lam,a,b): # lam in mm
    return brentq(lambda r: S(r,a)+S(r,b)-m*lam/D, 1e-9, 1e4)
b=150.0
# (2) paraxial
r_,a_,b_,lam_,m_=sp.symbols('r a b lambda m',positive=True)
print('(2)', sp.solve(sp.Eq(r_/a_+r_/b_, m_*lam_/d), r_))
# (3) a so that red(0.64um) first order at r=40
sb=S(40,b); sa=0.64e-3/D - sb; a=40*np.sqrt(1-sa**2)/sa
print('(3) sin b',sb,' sin a',sa,' a=',a)
print('   paraxial a:', 1/(0.64e-3/(D*40)-1/b))
# (4)
for m,lam in [(1,0.40e-3),(1,0.64e-3),(2,0.40e-3),(2,0.64e-3)]:
    print('(4) m',m,'lam',lam,' r=',ring(m,lam,a,b))
# (5) a -> infinity
print('(5) parallel red r=', b*np.tan(np.arcsin(0.4)), ' violet', b*np.tan(np.arcsin(0.25)))
print('    a=b=150 red', ring(1,0.64e-3,150,150))
# 追加：r=58 mm に来る波長
lam2=(S(58,a)+S(58,b))*D/2; print('second order lam at r=58 (um)',lam2*1e3)
lamp=S(58,b)*D; print('parallel first order lam at r=58 (um)',lamp*1e3)
# 数値で光路差を直接計算して確認：隣り合う溝 r, r+D からの S->溝->Q 経路差
r=40.0; L=lambda rr: np.hypot(rr,a)+np.hypot(rr,b)
print('path diff / lam (should be 1):',(L(r+D)-L(r))/0.64e-3)
