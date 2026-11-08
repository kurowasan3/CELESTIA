# Q200-1 成長しながら上昇する泡：v ∝ r^2, r = r0 + kt -> h ∝ r^3
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
rho,g,eta,k,r,r0,t,h,T=sp.symbols('rho g eta k r r0 t h T',positive=True)
v=2*rho*g*r**2/(9*eta)
H=sp.integrate(v/k,(r,r0,r)); print('h(r)=',sp.factor(H))
acc=sp.diff(v.subs(r,r0+k*t),t).subs(t,(r-r0)/k); print('a=',sp.simplify(acc))
rh=sp.solve(sp.Eq(H.subs(r0,0),h),r)[0]; print('r(h)=',rh); vh=sp.simplify(v.subs(r,rh)); print('v(h)=',vh)
vals={rho:1.0e3,g:9.8,eta:0.10,k:1.0e-5}
print('h at r=0.50mm (r0=0):',sp.N(H.subs(r0,0).subs(vals).subs(r,0.5e-3),3))
print('v at r=0.5mm',sp.N(v.subs(vals).subs(r,0.5e-3),3),' Re',sp.N((rho*v*2*r/eta).subs(vals).subs(r,0.5e-3),2))
print('a at r=0.5mm',sp.N(acc.subs(vals).subs(r,0.5e-3),3))
# 連なった泡の間隔 ∝ h^{2/3}: 1cm で 2.0mm → 8cm で
print('spacing at 8cm',2.0*8**(2/3))
# 数値：慣性を入れた運動方程式（付加質量 1/2 ρ V）でも準静的近似が成り立つか
pv=dict(rho=1.0e3,g=9.8,eta=0.10,k=1.0e-5)
def f(tt,s):
    rr=1e-6+pv['k']*tt; V=4/3*np.pi*rr**3; madd=0.5*pv['rho']*V
    # d(madd v)/dt = buoy - drag
    dmadd=0.5*pv['rho']*4*np.pi*rr**2*pv['k']
    return [s[1],(pv['rho']*V*pv['g']-6*np.pi*pv['eta']*rr*s[1]-dmadd*s[1])/madd]
sol=solve_ivp(f,[0,49.9],[0,0],rtol=1e-9,atol=1e-14,method='LSODA')
rr=1e-6+pv['k']*49.9; print('num v',sol.y[1,-1],' quasi',2*pv['rho']*pv['g']*rr**2/(9*pv['eta']),' num h',sol.y[0,-1],' quasi', 2*pv['rho']*pv['g']*(rr**3-1e-18)/(27*pv['eta']*pv['k']))
# 圧力効果：深さ 0.10 m を上がる間の半径変化（等温）
v1=float(vh.subs(vals).subs(h,0.010)); print('v(1cm)',v1,' T=',2.0e-3/v1, ' v(8cm)',float(vh.subs(vals).subs(h,0.08)))
print('r at 1cm',float(rh.subs(vals).subs(h,0.01)))

# --- (4) 気体を取り込まず等温膨張だけの場合：v ∝ r^2 ∝ P^{-2/3}
P0s,Ds,rhos,gs=sp.symbols('P0 D rho g',positive=True)
z=sp.symbols('z',positive=True)
r_of_z=(P0s+rhos*gs*z)**sp.Rational(-1,3)          # r ∝ P^{-1/3} (PV=const)
ratio=sp.simplify((r_of_z.subs(z,0)/r_of_z.subs(z,Ds))**2)
print('v_s/v_D =',ratio)
assert sp.simplify(ratio-(1+rhos*gs*Ds/P0s)**sp.Rational(2,3))==0
x=sp.symbols('x'); gam=sp.series((1+x)**sp.Rational(2,3),x,0,2).removeO().coeff(x,1); print('gamma',gam); assert gam==sp.Rational(2,3)
D2=sp.solve(sp.Eq((1+rhos*gs*Ds/P0s)**sp.Rational(2,3),2),Ds)[0]; print('D2 =',sp.simplify(D2))
assert sp.simplify(D2-(2*sp.sqrt(2)-1)*P0s/(rhos*gs))==0
print('D2 numeric (rho=1e3,P0=1e5,g=9.8):',float(D2.subs({P0s:1.0e5,rhos:1.0e3,gs:9.8})))
# 成長モデル：v ∝ h^{2/3} で速さ2倍 ⇔ h が 2^{3/2} 倍
hh=sp.symbols('hh',positive=True); fac=sp.solve(sp.Eq((hh)**sp.Rational(2,3),2),hh)[0]; print('height factor',fac,float(fac))
assert sp.simplify(fac-2*sp.sqrt(2))==0
# 参考：深さ 10 cm からの比
Dn=0.10; P0n=1.0e5
print('ratio for D=10cm:',float(ratio.subs({P0s:P0n,rhos:1.0e3,gs:9.8,Ds:Dn})))
print('(4) OK')
