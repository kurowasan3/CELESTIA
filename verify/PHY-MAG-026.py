# PHY-MAG-026 調和型中心力＋一様磁場：2つの円運動の角速度・回転系・静かに放した粒子の軌道
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,k,q,B,w=sp.symbols('m k q B omega',positive=True)
# clockwise (viewed from +z, B along +z, q>0): Lorentz inward
w1=sp.solve(sp.Eq(m*w**2,k+q*B*w),w); w2=sp.solve(sp.Eq(m*w**2,k-q*B*w),w)
w1=[r for r in w1 if r.subs({m:1,q:1,B:1,k:2})>0][0]; w2=[r for r in w2 if r.subs({m:1,q:1,B:1,k:2})>0][0]
print('cw',w1,' ccw',w2)
print('diff',sp.simplify(w1-w2),' prod',sp.simplify(w1*w2),' mean',sp.simplify((w1+w2)/2))
# numeric: m=q=B=1 -> wc=1 ; k=2 (=> w1=2, w2=1)
mm,qq,BB,kk=1.0,1.0,1.0,2.0
print('w1,w2 numeric', float(w1.subs({m:1,q:1,B:1,k:2})), float(w2.subs({m:1,q:1,B:1,k:2})))
def f(t,y):
    x,yy,vx,vy=y
    ax=(-kk*x+qq*(vy*BB))/mm; ay=(-kk*yy-qq*(vx*BB))/mm
    return [vx,vy,ax,ay]
r0=1.0
s=solve_ivp(f,[0,2*np.pi*3],[r0,0,0,0],rtol=1e-11,atol=1e-13,dense_output=True,max_step=0.001)
r=np.hypot(s.y[0],s.y[1]); sp_=np.hypot(s.y[2],s.y[3])
print('rmin num',r.min(),' theory 1/3')
i=np.argmin(r[s.t<2]); print('first rmin time',s.t[s.t<2][i],' theory pi/(2w\')=',np.pi/(2*1.5))
# stops (speed zero) in one period 2pi/wc
tt=s.t[(s.t>0.01)&(s.t<2*np.pi+0.01)]; vv=sp_[(s.t>0.01)&(s.t<2*np.pi+0.01)]
from scipy.signal import argrelmin
mins=[j for j in argrelmin(vv)[0] if vv[j]<1e-3]; print('stops in one period:',len(mins),'times',tt[mins])
yend=s.sol(2*np.pi); print('state after 2pi',yend)
# general r_min formula check with k=3
kk=3.0
s=solve_ivp(f,[0,40],[r0,0,0,0],rtol=1e-11,atol=1e-13,max_step=0.002)
wc=1; wp=np.sqrt(kk+0.25)
print('k=3 rmin num',np.hypot(s.y[0],s.y[1]).min(),' theory',0.5/wp)
