# K20151-1 verification
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
th,g,l,W=sp.symbols('theta g ell Omega',positive=True)
f=-g/l*sp.sin(th)+W**2*sp.sin(th)*sp.cos(th)   # theta'' in rotating frame
eq=sp.solve(sp.Eq(f,0),th); print('equilibria:',eq)
w2_low=sp.simplify(-sp.diff(f,th).subs(th,0)); print('omega^2 (theta0=0):',w2_low)
th0=sp.acos(g/(l*W**2)); w2_up=sp.simplify(-sp.diff(f,th).subs(th,th0)); print('omega^2 upper:',w2_up, ' = W^2 sin^2? ',sp.simplify(w2_up-W**2*sp.sin(th0)**2))
print('omega=Omega lower:',sp.solve(sp.Eq(w2_low,W**2),W))
print('omega=Omega upper:',sp.solve(sp.Eq(w2_up,W**2),W))
# numeric: g=9.8,l=1, Omega=sqrt(g/2l), a=0.01
gn,ln=9.8,1.0
def run(Wn,a=0.01,T=None):
    def rhs(t,y): return [y[1], -gn/ln*np.sin(y[0])+Wn**2*np.sin(y[0])*np.cos(y[0])]
    T=T or 40/ Wn
    s=solve_ivp(rhs,[0,T],[a,0],rtol=1e-11,atol=1e-13,dense_output=True,max_step=0.01)
    t=np.linspace(0,T,20001); th_,thd=s.sol(t)
    X=ln*np.sin(th_)*np.cos(Wn*t); Y=ln*np.sin(th_)*np.sin(Wn*t)
    vx=np.gradient(X,t); vy=np.gradient(Y,t)
    sp_=np.hypot(vx,vy)
    Fperp=2*1.0*Wn*ln*thd*np.cos(th_)   # m=1 ; d(m rho^2 Omega)/dt / rho
    return t,X,Y,sp_[5:-5],Fperp
Wn=np.sqrt(gn/(2*ln)); a=0.01
t,X,Y,spd,F=run(Wn,a)
r=np.hypot(X-a*ln/2,Y)
print('circle radius range:',r.min(),r.max(),' expected',a*ln/2)
print('speed range:',spd.min(),spd.max(),' expected aℓΩ=',a*ln*Wn)
print('max Fperp:',np.abs(F).max(),' expected m a g=',a*gn)
# general case Omega=0.5*Omega_c
Wn=0.5*np.sqrt(gn/ln); wn=np.sqrt(gn/ln-Wn**2)
t,X,Y,spd,F=run(Wn,a,T=200)
print('general: speed min/max',spd.min(),spd.max(),' expected',a*ln*min(Wn,wn),a*ln*max(Wn,wn))
