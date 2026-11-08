# J103-1: field lines of +4Q at A(0,0) and -Q at B(l,0); flux-tube conservation
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
l=1.0; QA=4.0; QB=-1.0
def E(p):
    x,y=p; r1=np.hypot(x,y); r2=np.hypot(x-l,y)
    return np.array([QA*x/r1**3+QB*(x-l)/r2**3, QA*y/r1**3+QB*y/r2**3])
def rhs(s,p):
    Ev=E(p); return Ev/np.linalg.norm(Ev)
def hitB(s,p): return np.hypot(p[0]-l,p[1])-1e-3
hitB.terminal=True
def far(s,p): return np.hypot(p[0],p[1])-400
far.terminal=True
def trace(alpha,r0=1e-3):
    p0=[r0*np.cos(alpha), r0*np.sin(alpha)]
    sol=solve_ivp(rhs,[0,2000],p0,events=[hitB,far],rtol=1e-10,atol=1e-12,max_step=0.5)
    return sol
# null point
xs=sp.symbols('x',positive=True); print('null:',sp.solve(sp.Eq(4/xs**2,1/(xs-1)**2),xs))
for adeg in [30,41.41,55,59,61,90,120,150]:
    a=np.radians(adeg); sol=trace(a)
    x,y=sol.y[0,-1],sol.y[1,-1]
    if sol.t_events[0].size:
        # arrival angle at B measured from direction B->A (i.e. -x): phi
        phi=np.degrees(np.arctan2(y, -(x-l)))
        pred=np.degrees(np.arccos(1-4*(1-np.cos(a))))
        print(f'alpha={adeg}: ends at B, phi={phi:.3f}, pred={pred:.3f}')
    else:
        # asymptotic direction measured from center of charge (-l/3,0)
        g=np.degrees(np.arctan2(y, x+l/3)); g2=np.degrees(np.arctan2(y,x))
        c=1-(4*(1-np.cos(a))-2)/3; pred=np.degrees(np.arccos(c))
        print(f'alpha={adeg}: to infinity, gamma(from centroid)={g:.3f} (from A {g2:.3f}), pred={pred:.3f}')
print('cos alpha=3/4 ->',np.degrees(np.arccos(0.75)))
