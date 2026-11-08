# M44-1: circular orbit r0=3R; probe (m1) ejected radially outward; shuttle (m-m1) gets inward radial speed v0/2.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
G,M,R,m,m1,r=sp.symbols('G M R m m1 r',positive=True)
r0=3*R; v0=sp.sqrt(G*M/r0); T0=2*sp.pi*r0/v0
print('(1) v0',v0,' T0',sp.simplify(T0))
w2=v0/2
h=r0*v0  # tangential component unchanged
E=sp.Rational(1,2)*(v0**2+w2**2)-G*M/r0
rs=sp.solve(sp.Eq(h**2/(2*r**2)-G*M/r,E),r); print('(2) apsides',rs)
vp=h/min(rs,key=lambda e: e.subs(R,1)); print('    v_perigee',sp.simplify(vp))
a=sum(rs)/2; print('(3) T/T0=',sp.simplify((a/r0)**sp.Rational(3,2)))
w1=(m-m1)*w2/m1
cond=sp.solve(sp.Eq(v0**2+w1**2,2*G*M/r0),m1); print('(4) boundary m1=',cond)
vinf=sp.sqrt(v0**2+w1**2-2*G*M/r0).subs(m1,m/4); print('(5) vinf',sp.simplify(vinf))
# numeric orbit (GM=1,R=1)
def f(t,y):
    x,yy,vx,vy=y; rr=np.hypot(x,yy); return [vx,vy,-x/rr**3,-yy/rr**3]
r0n=3.0; v0n=np.sqrt(1/r0n)
sol=solve_ivp(f,[0,200],[r0n,0,-v0n/2,v0n],rtol=1e-12,atol=1e-12,max_step=0.01)
rr=np.hypot(sol.y[0],sol.y[1]); print('num rmin,rmax',rr.min(),rr.max())
# period: angle wraps
ang=np.unwrap(np.arctan2(sol.y[1],sol.y[0])); idx=np.searchsorted(ang,2*np.pi)
T0n=2*np.pi*r0n/v0n
print('num T/T0', sol.t[idx]/T0n, ' th', float((sp.Rational(4,3))**sp.Rational(3,2)))
