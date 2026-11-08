# K20081-1: べき乗の中心力による円運動に小さな一定の抵抗力がはたらくとき
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
k,r,m,D,n=sp.symbols('k r m D n',real=True)
K=k*r**(n+1)/2            # from m v^2/r = k r^n
dVdr=k*r**n               # attractive: dV/dr = F
dKdr=sp.diff(K,r); dEdr=dKdr+dVdr
ratio=sp.simplify(dKdr/dEdr); print('dK/dE =',ratio)
# speed change rate: m v dv/dt = ratio * (-D v)
print('dv/dt =',sp.simplify(-ratio*D/m))
# numeric
def sim(nn,kk=1.0,mm=1.0,DD=1e-3,r0=1.0,tend=None):
    v0=np.sqrt(kk*r0**(nn+1)/mm)
    def rhs(t,y):
        x,yv,vx,vy=y; rr=np.hypot(x,yv); sp_=np.hypot(vx,vy)
        F=kk*rr**nn
        return [vx,vy,-F*x/rr/mm-DD*vx/sp_/mm,-F*yv/rr/mm-DD*vy/sp_/mm]
    s=solve_ivp(rhs,[0,tend],[r0,0,0,v0],rtol=1e-10,atol=1e-12,dense_output=True,max_step=0.05)
    return s,v0
for nn in [1.0,-2.0,-1.0,-1.5,0.0]:
    tend=60.0
    s,v0=sim(nn,tend=tend)
    ts=np.linspace(0,tend,6001); Y=s.sol(ts); spd=np.hypot(Y[2],Y[3]); rr=np.hypot(Y[0],Y[1])
    slope=np.polyfit(ts,spd,1)[0]
    th=-(nn+1)/(nn+3)*1e-3
    print('n=%5.1f  numeric dv/dt=%.4e  theory=%.4e   r: %.4f -> %.4f'%(nn,slope,th,rr[0],rr[-1]))
# n=-1: r(t)=r0 exp(-D t/(m v)), v=sqrt(k/m)
s,v0=sim(-1.0,tend=200.0); ts=np.linspace(0,200,5); Y=s.sol(ts); rr=np.hypot(Y[0],Y[1])
print('n=-1 r num',rr,' th',np.exp(-1e-3*ts/(1.0*v0)))
# gravity: v grows D/m; r=GM/v^2 ; time for r -> r0/4 is m v0/D
s,v0=sim(-2.0,tend=1000.0); ts=np.array([0,500,1000.]); Y=s.sol(ts)
print('n=-2 speed num',np.hypot(Y[2],Y[3]),' th',v0+1e-3*ts)
