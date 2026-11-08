# Q075-1: 大きな水槽につないだ水平管の流れの立ち上がり
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,H,L,t,rho=sp.symbols('g H L t rho',positive=True)
V=sp.sqrt(2*g*H); tau=2*L/V
v=V*sp.tanh(t/tau)
# エネルギー収支 L dv/dt = gH - v^2/2
assert sp.simplify(L*sp.diff(v,t)-(g*H-v**2/2))==0
print('tau =',sp.simplify(tau))
# 数値
gv,Hv,Lv=9.8,5.0,20.0
Vv=np.sqrt(2*gv*Hv); tv=2*Lv/Vv
s=solve_ivp(lambda tt,y:[(gv*Hv-y[0]**2/2)/Lv, y[0]],[0,20],[0,0],rtol=1e-11,atol=1e-12,dense_output=True)
t90=tv*np.arctanh(0.9)
print('V=%.4f tau=%.4f t90=%.4f'%(Vv,tv,t90), 'v(t90)/V=',s.sol(t90)[0]/Vv)
assert abs(s.sol(t90)[0]/Vv-0.9)<1e-8
# t90 までに出口から出た水柱の長さ ∫v dt = 2L ln cosh(t/tau)
X=2*Lv*np.log(np.cosh(t90/tv)); print('X=%.4f num=%.4f'%(X,s.sol(t90)[1])); assert abs(X-s.sol(t90)[1])<1e-6
print('2L ln(1/sqrt(1-0.81)) =',2*Lv*np.log(1/np.sqrt(0.19)))
# 管の入口の圧力（ゲージ圧）: rho L dv/dt = rho(gH - v^2/2)
x=sp.symbols('x',positive=True)
pin=rho*L*sp.diff(v,t)
print('p_in - p0 at v=V/2:', sp.simplify((rho*(g*H-(V/2)**2/2))))
# 半分の速さになる時刻から L を逆算: tanh(t/tau)=1/2 -> t=tau*atanh(1/2)
print('atanh(1/2)=',np.arctanh(0.5),' L for t_half=3.0 s, H=5:', 3.0*Vv/(2*np.arctanh(0.5)))
print('OK')
