# K20123-1 verification: clock rates (SR + gravitational potential), first order in 1/c^2
import sympy as sp
G,M,R,r,c,W,eps,Om=sp.symbols('G M R r c W epsilon Omega',positive=True)
# (1) rotating disk: inertial SR vs centrifugal potential
rp=sp.symbols('rp',positive=True)
sr = 1 - (W*rp)**2/(2*c**2)
phi = lambda x: -sp.Rational(1,2)*W**2*x**2
gr = 1 + (phi(rp)-phi(0))/c**2
print('(1)', sp.simplify(sr-gr))
# (2) satellite vs ground (inertial frame, no rotation): rate factor relative to clock at infinity
rate = lambda pot,u2: 1 + pot/c**2 - u2/(2*c**2)
sat = rate(-G*M/r, G*M/r)
gnd = rate(-G*M/R, 0)
ratio = sp.series(sat/gnd, c, sp.oo, 3).removeO()
ratio = sp.expand(sat-gnd)  # first order equivalent
print('(2) ratio-1 =', sp.simplify(ratio))
print('(3) r0 =', sp.solve(sp.Eq(ratio,0), r))
# (4) with Earth rotation
gnd2 = rate(-G*M/R, (Om*R)**2)
r0p = sp.solve(sp.Eq(sat-gnd2,0), r)[0]
print('(4) r0 exact-in-first-order =', sp.simplify(r0p))
d=sp.symbols('delta',positive=True)
expr = sp.simplify(r0p.subs(Om, sp.sqrt(d*G*M/R**3)))
print('    in delta:', expr, ' series:', sp.series(expr, d, 0, 2))
g=9.8; Rn=6.4e6; Omn=7.3e-5; cn=3.0e8
dn = Omn**2*Rn/g
print('    delta=',dn,' shift km =', 1.5*Rn*dn/2/1e3, ' exact shift km=', (1.5*Rn - 3*Rn/(2*(1+dn/2)))/1e3)
# ---- (5) elliptical orbit, semi-major axis a (no Earth rotation) ----
a=sp.symbols('a',positive=True)
v2=G*M*(2/r-1/a)
satE=rate(-G*M/r, v2)
ratioE=sp.expand(satE-gnd)   # first order
print('(5) ratio-1 =',sp.simplify(ratioE))
assert sp.simplify(ratioE-G*M/c**2*(1/R-2/r+1/(2*a)))==0
assert sp.simplify(ratioE.subs(a,r)-sp.expand(sat-gnd))==0
# numeric: integrate Kepler orbit, check <1/r>=1/a and mean rate equals circular radius a
import numpy as np
from scipy.integrate import solve_ivp
GMv=1.0; rp_=1.0; e_=0.6; av=rp_/(1-e_); vp=np.sqrt(GMv*(1+e_)/rp_)
def rhs(t,y):
    x,yy,vx,vy=y; d=np.hypot(x,yy); return [vx,vy,-GMv*x/d**3,-GMv*yy/d**3]
T=2*np.pi*np.sqrt(av**3/GMv)
so=solve_ivp(rhs,[0,T],[rp_,0,0,vp],rtol=1e-12,atol=1e-12,dense_output=True)
ts=np.linspace(0,T,200001); Y=so.sol(ts); rr=np.hypot(Y[0],Y[1]); vv2=Y[2]**2+Y[3]**2
inv=np.trapezoid(1/rr,ts)/T
print('(5) <1/r>=%.8f 1/a=%.8f'%(inv,1/av)); assert abs(inv-1/av)<1e-6
Rv=0.5; cc=1.0  # units: rates measured as coefficient of 1/c^2
inst=(1/Rv) - 1/rr - vv2/2      # (ratio-1)*c^2/GM with GM=1 (A at rest at R)
print('(5) mean rate coeff=%.8f  circular a: %.8f'%(np.trapezoid(inst,ts)/T, 1/Rv-1.5/av))
assert abs(np.trapezoid(inst,ts)/T-(1/Rv-1.5/av))<1e-6
imax=np.argmax(inst); print('(5) fastest at r=%.4f (apogee %.4f)'%(rr[imax],av*(1+e_)))
assert abs(rr[imax]-av*(1+e_))<1e-3
