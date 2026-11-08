# P17 類題：太陽に落とすのは脱出より難しい／いったん遠ざかってから落とす方が安い
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
w,V,rE,rp,RS=sp.symbols('w V r_E r_p R_S',positive=True)
# r_p/r_E from aphelion speed w at r_E (GM = V^2 rE)
vp=sp.symbols('vp',positive=True)
sol=sp.solve([sp.Eq(w*rE,vp*rp),sp.Eq(w**2/2-V**2,vp**2/2-V**2*rE/rp)],[vp,rp],dict=True)
for s_ in sol: print(sp.simplify(s_[rp]/rE))
VE=30.0; ve=11.2; ratio=1/215
u_esc=(np.sqrt(2)-1)*VE; v0_esc=np.hypot(ve,u_esc)
print('(1) u=',u_esc,' v0=',v0_esc)
wmax=VE*np.sqrt(2*ratio/(1+ratio)); u_dir=VE-wmax; v0_dir=np.hypot(ve,u_dir)
print('(2) wmax=',wmax,' u=',u_dir,' v0=',v0_dir)
# (3) 遠日点 8 rE の楕円
k=8
w1=VE*np.sqrt(2*k/(1+k)); u1=w1-VE; v01=np.hypot(ve,u1); wa=w1/k
Va=VE/np.sqrt(k); wamax=Va*np.sqrt(2*ratio/(k+ratio)); dv2=wa-wamax
print('(3) w1=',w1,' u1=',u1,' v0=',v01,' wa=',wa,' wamax=',wamax,' dv2=',dv2,' total=',v01+dv2)
# 数値軌道（単位 rE=1, VE=1, GM=1）
def orbit(vt,r0=1.0,T=200):
    f=lambda t,y:[y[2],y[3],-y[0]/np.hypot(y[0],y[1])**3,-y[1]/np.hypot(y[0],y[1])**3]
    s=solve_ivp(f,[0,T],[r0,0,0,vt],rtol=1e-11,atol=1e-13,max_step=0.01)
    r=np.hypot(s.y[0],s.y[1]); return r.min(),r.max()
print('direct: perihelion with w=wmax:',orbit(wmax/VE,T=4)[0],' target',ratio)
print('transfer: aphelion with w1:',orbit(w1/VE,T=40)[1])
print('from 8rE with wamax: perihelion',orbit(wamax/VE,r0=8.0,T=60)[0],' target',ratio)
assert abs(orbit(wmax/VE,T=4)[0]-ratio)<1e-5
