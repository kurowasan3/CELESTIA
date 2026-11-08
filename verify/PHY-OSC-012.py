# M40-1: bead in straight smooth tube, tube makes angle theta with vertical axis, rotates at omega.
# spring from O (on axis, lower end) natural length l, const k. r = distance along tube.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,k,g,l,w,th,t=sp.symbols('m k g l omega theta t',positive=True)
r=sp.Function('r')(t)
# Lagrangian in inertial frame (independent of centrifugal-force reasoning)
Tk=sp.Rational(1,2)*m*(r.diff(t)**2+w**2*r**2*sp.sin(th)**2)
V=sp.Rational(1,2)*k*(r-l)**2+m*g*r*sp.cos(th)
Lg=Tk-V
eom=sp.simplify(sp.diff(Lg.diff(r.diff(t)),t)-Lg.diff(r))
print('EOM: m r\'\' =', sp.simplify(-(eom-m*r.diff(t,2))))
K=k-m*w**2*sp.sin(th)**2
r0=(k*l-m*g*sp.cos(th))/K
print('(1) r0 check', sp.simplify((-(eom-m*r.diff(t,2))).subs(r,r0)))
T0=2*sp.pi*sp.sqrt(m/k); Tw=2*sp.pi*sp.sqrt(m/K)
# (3) release at r=l (natural length) from rest rel. to tube: amplitude |r0-l|
A=sp.simplify(r0-l); print('(3) r0-l =',sp.factor(A))
rmax=sp.simplify(2*r0-l); print('    farthest r =',sp.factor(rmax))
# (4) Tw=2T0
wsol=sp.solve(sp.Eq(Tw**2,4*T0**2),w); print('(4) omega =',wsol)
ws=wsol[0]; r0s=sp.simplify(r0.subs(w,ws)); print('    r0 =',r0s,' rmax=',sp.simplify(2*r0s-l))
# numeric check of (3),(4)
sub={m:0.1,k:20.0,g:9.8,l:0.3,th:sp.pi/3}
wn=float(ws.subs(sub))
def f(tt,y): return [y[1], (-(sub[k])*(y[0]-sub[l]) + sub[m]*wn**2*y[0]*np.sin(np.pi/3)**2 - sub[m]*sub[g]*np.cos(np.pi/3))/sub[m]]
rr=solve_ivp(f,[0,3],[sub[l],0],rtol=1e-11,atol=1e-13,max_step=1e-4)
print('num rmax',rr.y[0].max(),' th',float((2*r0s-l).subs(sub)))
# period from zero crossings of velocity
v=rr.y[1]; tt=rr.t; idx=np.where((v[:-1]>0)&(v[1:]<=0))[0]
print('num period',np.diff(tt[idx])[:2],' th',float((Tw.subs(w,ws)).subs(sub)), ' T0',float(T0.subs(sub)))
