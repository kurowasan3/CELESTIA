# N55-1 betatron with core field: orbit radius R=a sqrt(beta-1), energy bookkeeping
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
a,beta,Bo,r,e,m,Bf,tau,t=sp.symbols('a beta B_o r e m B_f tau t',positive=True)
Phi=sp.pi*a**2*beta*Bo+sp.pi*(r**2-a**2)*Bo
R=sp.solve(sp.Eq(Phi,2*sp.pi*r**2*Bo),r); print('R=',R)
Rv=a*sp.sqrt(beta-1)
# energy check with linear ramp Bo=Bf t/tau
emf=sp.diff(Phi.subs(r,Rv).subs(Bo,Bf*t/tau),t); print('emf per turn',sp.simplify(emf))
N=sp.integrate(e*Bf*t/tau/(2*sp.pi*m),(t,0,tau)); print('N=',N)
K=(e*Bf*Rv)**2/(2*m); print('N*e*emf - K =',sp.simplify(N*e*emf-K))
# numeric: full 2D integration of electron in field Bz(r,t)=f(t)*(beta*B0 if r<a else B0) with induced E_phi(r)=-(1/2pi r) dPhi/dt
an,betan,Bfn,taun=0.05,5.0,0.02,2e-6
me,ee=9.11e-31,1.602e-19
Rn=an*np.sqrt(betan-1)
def Bz(rr,tt): return Bfn*tt/taun*(betan if rr<an else 1.0)
def dPhi(rr): # dPhi/dt within radius rr (rr>a)
    return np.pi*(an**2*betan+(rr**2-an**2))*Bfn/taun
def f(tt,s):
    x,y,vx,vy=s; rr=np.hypot(x,y)
    Ephi=-dPhi(rr)/(2*np.pi*rr)
    Ex,Ey=-Ephi*y/rr, Ephi*x/rr
    b=Bz(rr,tt)
    ax=(-ee/me)*(Ex+vy*b); ay=(-ee/me)*(Ey-vx*b)
    return [vx,vy,ax,ay]
t0=taun*0.1; v0=ee*Bz(Rn,t0)*Rn/me
# electron with B along +z rotates counterclockwise -> start at (R,0) moving +y
sol=solve_ivp(f,[t0,taun],[Rn,0,0,v0],rtol=1e-10,atol=1e-14,max_step=1e-10)
rr=np.hypot(sol.y[0],sol.y[1]); print('R',Rn,'r min/max',rr.min(),rr.max())
vend=np.hypot(sol.y[2,-1],sol.y[3,-1]); print('v end',vend,'expected',ee*Bfn*Rn/me)
