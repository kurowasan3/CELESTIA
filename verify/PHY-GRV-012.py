# O20211-1 : free tether co-rotating with the Earth; tension maximum; release from the tether
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
r,R0,Rs,lam,w,G,M,rp,rr=sp.symbols('r R_0 R_s lambda omega G M r_p rr',positive=True)
GM=w**2*Rs**3
# element: lam dr * w^2 r = lam dr GM/r^2 - dT  -> dT/dr = lam(GM/r^2 - w^2 r), T(R0)=0
T=sp.integrate(lam*(GM/rr**2-w**2*rr),(rr,R0,r))
print('T(r)=',sp.simplify(T))
print('dT/dr=0 at',sp.solve(sp.diff(T,r),r))
Tmax=sp.factor(sp.simplify(T.subs(r,Rs))); print('Tmax=',Tmax)
# top end: T(R2)=0
R2=sp.symbols('R_2',positive=True)
print('R2 eq:',sp.factor(sp.simplify(T.subs(r,R2))))
# release at r with speed w r (tangential). other apsis from L & E conservation
v=w*r
sol=sp.solve([sp.Eq(r*v,rp*sp.Symbol('vp')), sp.Eq(v**2/2-GM/r, sp.Symbol('vp')**2/2-GM/rp)],[rp,sp.Symbol('vp')],dict=True)
for s_ in sol: print('other apsis:',sp.simplify(s_[rp]))
s=sp.symbols('s',positive=True)
print('in s=(r/Rs)^3: r*s/(2-s) check:',sp.simplify(r*(r/Rs)**3/(2-(r/Rs)**3)-[q[rp] for q in sol if sp.simplify(q[rp]-r)!=0][0]))
# numeric: Rs=6.6 R0, minimal release radius below Rs that misses Earth
k=1/6.6
rho=brentq(lambda p:p**4+k*p**3-2*k,0.3,1.0); print('rho=r/Rs=',rho,' r/R0=',rho*6.6)
# orbit integration (units R0=1, GM=1)
GMv=1.0; Rsv=6.6; wv=np.sqrt(GMv/Rsv**3)
def peri(r0):
    def f(t,y):
        x,yy,vx,vy=y; d=np.hypot(x,yy)**3; return [vx,vy,-GMv*x/d,-GMv*yy/d]
    s=solve_ivp(f,[0,2e3],[r0,0,0,wv*r0],rtol=1e-11,atol=1e-12,max_step=0.5)
    return np.min(np.hypot(s.y[0],s.y[1])), np.max(np.hypot(s.y[0],s.y[1]))
for r0 in [rho*6.6,4.0,5.5,6.6,7.5]:
    s3=(r0/Rsv)**3; print('r0=%.3f numeric min/max r'%r0,peri(r0),' formula other apsis',r0*s3/(2-s3))
print('escape radius 2^(1/3)Rs =',2**(1/3)*6.6)
