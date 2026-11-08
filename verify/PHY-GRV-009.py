# K20191-1: phasing rendezvous
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
th,nn,eps=sp.symbols('theta n epsilon',positive=True)
r=1-th/(2*sp.pi*nn)          # T1/T0
a1=r**sp.Rational(2,3)       # a/R
d=2*a1-1
V=sp.sqrt(2-1/a1)
print('V1/V0 alt form check', sp.simplify(V-sp.sqrt(2*d/(1+d))))
# n=1, theta=pi
for n_ in [1,2]:
    print('decel n=',n_,'d/R=',float(d.subs({th:sp.pi,nn:n_})),'V1/V0=',float(V.subs({th:sp.pi,nn:n_})))
# min n with d>=1/2
nmin=next(n_ for n_ in range(1,10) if float(d.subs({th:sp.pi,nn:n_}))>=0.5); print('nmin',nmin)
ra=1+sp.Rational(1,4); a2=ra**sp.Rational(2,3); V2=sp.sqrt(2-1/a2)
print('accel V/V0',float(V2),'apo',float(2*a2-1))
dv_dec=2*(1-float(V.subs({th:sp.pi,nn:2}))); dv_acc=2*(float(V2)-1)
print('total dV dec',dv_dec,'acc',dv_acc)
f=sp.sqrt(2-(1+eps)**sp.Rational(-2,3)); print('series',sp.series(f,eps,0,3))
# numerical orbit check: GM=1,R=1
def rhs(t,y): x,yv,vx,vy=y; rr=(x*x+yv*yv)**1.5; return [vx,vy,-x/rr,-yv/rr]
T0=2*np.pi
for Vfac,nrev,T1 in [(float(V.subs({th:sp.pi,nn:2})),2,0.75*T0),(float(V2),2,1.25*T0)]:
    s=solve_ivp(rhs,[0,nrev*T1],[1,0,0,Vfac],rtol=1e-11,atol=1e-12,dense_output=True)
    ts=np.linspace(0,nrev*T1,20001); rs=np.hypot(*s.sol(ts)[:2])
    end=s.y[:,-1]
    # station starts ahead by pi: angle pi + t (omega=1)
    ang=np.pi+nrev*T1; zx,zy=np.cos(ang),np.sin(ang)
    print('end pos',end[:2],'station',(zx,zy),'rmin',rs.min(),'rmax',rs.max())
