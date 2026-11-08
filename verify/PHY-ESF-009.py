# K20232-1: hemispherical analyzer as Kepler problem
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
a,b,V0,e,m,r,r0,K0,dl,s=sp.symbols('a b V0 e m r r0 K0 delta s',positive=True)
# field between concentric spheres with potential difference V0: E = V0*a*b/((b-a) r^2)
Q=sp.symbols('Q',positive=True); k=sp.symbols('k',positive=True)
Vdiff=k*Q/a-k*Q/b
Qs=sp.solve(sp.Eq(Vdiff,V0),Q)[0]
E=sp.simplify(k*Qs/r**2); print('E=',E)
kap=e*V0*a*b/(b-a)     # force = kap/r^2
K0s=sp.simplify(kap/(2*r0)); print('K0=',K0s)
# exit radius for K=K0(1+delta), perpendicular entry
r1=sp.symbols('r1',positive=True)
v=sp.sqrt(2*K0s*(1+dl)/m); v1=v*r0/r1
eq=sp.Eq(m*v**2/2-kap/r0, m*v1**2/2-kap/r1)
sols=sp.solve(eq,r1); print('r1 sols',[sp.simplify(x) for x in sols])
# window for slit half-width w/2 = s*r0
d1=sp.solve(sp.Eq((1+dl)/(1-dl),1+s),dl); d2=sp.solve(sp.Eq((1+dl)/(1-dl),1-s),dl)
print('delta+',d1,' delta-',d2, ' width',sp.simplify(d1[0]-d2[0]))
# numeric orbit check
an,bn,r0n=1.0,1.6,1.3; kapn=1.0; mn=1.0
K0n=kapn/(2*r0n)
for dd in [0.05,-0.05,0.1]:
    vn=np.sqrt(2*K0n*(1+dd)/mn)
    def rhs(t,y):
        x,yy,vx,vy=y; rr=np.hypot(x,yy); A=-kapn/mn/rr**3
        return [vx,vy,A*x,A*yy]
    ev=lambda t,y: y[1]; ev.terminal=True; ev.direction=-1  # cross y=0 going to negative? start at (-r0,0) moving +y
    sol=solve_ivp(rhs,[0,100],[-r0n,0,0,vn],events=ev,rtol=1e-12,atol=1e-12)
    # start at X=(-r0,0) moving +y, clockwise around to Y on +x axis
    ye=sol.y_events[0][0]
    print(dd,'numeric exit r',ye[0],' theory',r0n*(1+dd)/(1-dd), ' max r', np.max(np.hypot(sol.y[0],sol.y[1])))
# photoemission numeric example: thresholds
