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

# ---- (3) general apoapsis k r_E ----
k_,eps_=sp.symbols('k epsilon',positive=True)
w1s,wq=sp.symbols('w1 wq',positive=True)
sol3=sp.solve([sp.Eq(w1s,k_*wq),sp.Eq(w1s**2/2-V**2,wq**2/2-V**2/k_)],[w1s,wq],dict=True)
sol3=[q for q in sol3 if sp.simplify(q[w1s]-V)!=0]
print('(3) w1 =',[sp.simplify(q[w1s]) for q in sol3])
w1f=V*sp.sqrt(2*k_/(k_+1)); wqf=V*sp.sqrt(2/(k_*(k_+1)))
assert any(sp.simplify(q[w1s]-w1f)==0 and sp.simplify(q[wq]-wqf)==0 for q in sol3)
# perihelion from Q with speed w' : rp/(k rE) = w'^2/(2VQ^2-w'^2), VQ^2=V^2/k ; set rp = eps rE
wp=sp.symbols('wp',positive=True)
wpmax=sp.solve(sp.Eq(wp**2/(2*V**2/k_-wp**2),eps_/k_),wp)[0]
assert sp.simplify(wpmax-V*sp.sqrt(2*eps_/(k_*(k_+eps_))))==0
dwq=wqf-wpmax
print('(3) k->oo: u ->',sp.limit(w1f-V,k_,sp.oo),'  dwQ ->',sp.limit(dwq,k_,sp.oo))
ve_,VE_=11.2,30.0
tot=lambda kk: np.hypot(ve_,VE_*(np.sqrt(2*kk/(kk+1))-1))+VE_*(np.sqrt(2/(kk*(kk+1)))-np.sqrt(2*ratio/(kk*(kk+ratio))))
ks=np.exp(np.linspace(np.log(1.01),np.log(1e6),400)); tv=tot(ks)
assert np.all(np.diff(tv)<0); print('(3) total decreasing in k; k=8:',tot(8.0),' k=1e6:',tot(1e6),' escape v0:',v0_esc)
assert abs(tot(8.0)-(v01+dv2))<1e-9
# numeric orbit: from 1 rE with w1(k) reaches k ; from k with wpmax reaches eps (units rE=1, VE=1)
for kk in [3.0,20.0]:
    w1n=np.sqrt(2*kk/(kk+1)); wpn=np.sqrt(2*ratio/(kk*(kk+ratio)))
    print('(3) k=%g aphelion %.6f ; perihelion from Q %.6f (target %.6f)'%(kk,orbit(w1n,T=4*np.pi*((kk+1)/2)**1.5)[1],orbit(wpn,r0=kk,T=2*np.pi*(kk/2)**1.5*1.2)[0],ratio))
# time to Q: half period with a=(k+1)/2 rE
print('(3) time to Q for k=8: %.3f T_E'%(0.5*((8+1)/2)**1.5))
