# Q075-1: 大きな水槽につないだ水平管の流れの立ち上がり
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,H,L,t,rho=sp.symbols('g H L t rho',positive=True)
V=sp.sqrt(2*g*H); tau=2*L/V
v=V*sp.tanh(t/tau)
# エネルギー収支 L dv/dt = gH - v^2/2
assert sp.simplify(L*sp.diff(v,t)-(g*H-v**2/2))==0
print('tau =',sp.simplify(tau))
# 管の入口の圧力（ゲージ圧）: rho L dv/dt = rho(gH - v^2/2)
pin=rho*L*sp.diff(v,t)
assert sp.simplify(pin-rho*(g*H-v**2/2))==0
print('p_in - p0 = rho(gH - v^2/2) OK')

# ---- (5) 太さの変わる管: L1,S1 (水槽側) + L2,S2 (出口側) ----
L1,L2,S1,S2,vv,a=sp.symbols('L1 L2 S1 S2 v a',positive=True)
v1=S2/S1*vv
# エネルギー収支: rho S2 v gH = d/dt(KE_tube) + 1/2 rho S2 v^3, with dv/dt = a
KE=sp.Rational(1,2)*rho*S1*L1*v1**2+sp.Rational(1,2)*rho*S2*L2*vv**2
dKE=sp.diff(KE,vv)*a
a_energy=sp.solve(sp.Eq(rho*S2*vv*g*H, dKE+sp.Rational(1,2)*rho*S2*vv**3),a)[0]
Leff=L2+S2/S1*L1
assert sp.simplify(a_energy-(g*H-vv**2/2)/Leff)==0
# 運動方程式ルート: 入口 p0+rho g H - rho v1^2/2, 変わり目 p1 + rho v1^2/2 = p2 + rho v^2/2
p0,p1=sp.symbols('p0 p1')
p2=p1+rho*v1**2/2-rho*vv**2/2
pin2=p0+rho*g*H-rho*v1**2/2
# rho L1 dv1/dt = pin2 - p1 ; rho L2 dv/dt = p2 - p0 ; dv1/dt = S2/S1 a
sol=sp.solve([sp.Eq(rho*L1*S2/S1*a, pin2-p1), sp.Eq(rho*L2*a, p2-p0)],[a,p1],dict=True)[0]
assert sp.simplify(sol[a]-a_energy)==0
print('(5) energy and momentum routes agree: Leff =',Leff)
# 終端速度と時定数
tau2=2*Leff/V
vt=V*sp.tanh(t/tau2)
assert sp.simplify(Leff*sp.diff(vt,t)-(g*H-vt**2/2))==0
print('(5) V\' = sqrt(2gH), tau\' =',sp.simplify(tau2))
# 極限
assert sp.simplify(Leff.subs(S1,S2)-(L1+L2))==0
assert sp.limit(Leff,S1,sp.oo)==L2
# 数値: ODE を直接積分して tanh 解と比較
gv,Hv,L1v,L2v,S1v,S2v=9.8,5.0,12.0,8.0,3.0,1.0
Le=L2v+S2v/S1v*L1v; Vv=np.sqrt(2*gv*Hv); tv=2*Le/Vv
s=solve_ivp(lambda tt,y:[(gv*Hv-y[0]**2/2)/Le],[0,30],[0.0],rtol=1e-11,atol=1e-12,dense_output=True)
for tt in [1.0,3.0,10.0]:
    assert abs(s.sol(tt)[0]-Vv*np.tanh(tt/tv))<1e-7
print('(5) numeric check OK, Leff=%.3f tau=%.3f'%(Le,tv))
print('OK')
