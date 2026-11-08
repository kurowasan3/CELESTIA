# Q064-1: caterpillar on uniformly stretched string: (a) winch reels in at fixed position, (b) end moves away
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
t,L,v0,c,x0=sp.symbols('t L v0 c x0',positive=True)
X=sp.Function('X')
# case (a): dx/dt = v0 x/L - c
sol=sp.dsolve(sp.Eq(X(t).diff(t),v0*X(t)/L-c),X(t),ics={X(0):x0}); print(sol)
ta=sp.solve(sp.Eq(sol.rhs.subs(x0,c*L/(2*v0)),0),t); print('time (a), x0=cL/2v0:',ta)
# case (b): L(t)=L+v0 t, dx/dt = v0 x/(L+v0 t) - c
solb=sp.dsolve(sp.Eq(X(t).diff(t),v0*X(t)/(L+v0*t)-c),X(t),ics={X(0):x0}); print(sp.simplify(solb.rhs))
tb=sp.solve(sp.Eq(sp.simplify(solb.rhs),0),t); print('time (b):',tb)
# numeric sanity of (3) formula (L=1, v0=0.01, c=0.001, x0=0.5)
Ln,vn,cn,xn=1.0,0.01,0.001,0.5
ev=lambda tt,y:y[0]; ev.terminal=True
r=solve_ivp(lambda tt,y:[vn*y[0]/(Ln+vn*tt)-cn],[0,1e6],[xn],events=ev,rtol=1e-10,atol=1e-12)
print('(3) numeric time',r.t_events[0], ' formula',(Ln/vn)*(np.exp(vn*xn/(cn*Ln))-1))
assert abs(r.t_events[0][0]-(Ln/vn)*(np.exp(vn*xn/(cn*Ln))-1))<1e-3*r.t_events[0][0]
# ---- (4) l(t)=L(1+t/tau)^n ----
tau,n,ss=sp.symbols('tau n s',positive=True)
ell=L*(1+t/tau)**n
u=X(t)/ell
du=sp.simplify(sp.diff(u,t).subs(sp.Derivative(X(t),t),sp.diff(ell,t)/ell*X(t)-c))
print('(4) d(x/l)/dt =',du); assert sp.simplify(du+c/ell)==0
Tn=tau*((1+(1-n)*x0/(c*tau))**(1/(1-n))-1)
# check: integral_0^T c/l dt = x0/L
for nv in [sp.Rational(1,2),sp.Rational(1,3),sp.Rational(3,4)]:
    I_=sp.integrate((c/ell).subs(n,nv),(t,0,Tn.subs(n,nv)))
    assert sp.simplify(I_-x0/L)==0
print('(4) T(n<1) satisfies integral condition')
print('(4) n->1 limit:',sp.simplify(sp.limit(Tn,n,1)))
assert sp.simplify(sp.limit(Tn,n,1)-tau*(sp.exp(x0/(c*tau))-1))==0
print('(4) n>1 total integral:',sp.simplify(sp.integrate((c/ell).subs(n,2),(t,0,sp.oo))),' general c tau/(L(n-1))')
# numeric: n=0.5 reach time, n=2 threshold
def run(nv,tauv,x0v,tmax):
    f=lambda tt,y:[ (nv/tauv)/(1+tt/tauv)*y[0]-cn ]
    rr=solve_ivp(f,[0,tmax],[x0v],events=ev,rtol=1e-11,atol=1e-13)
    return rr.t_events[0]
tauv=100.0
tn=run(0.5,tauv,0.5,1e7); th=tauv*((1+0.5*0.5/(cn*tauv))**2-1)
print('(4) n=0.5 numeric T=%.3f formula %.3f'%(tn[0],th)); assert abs(tn[0]-th)<1e-4*th
thr=cn*tauv/(2-1)
print('(4) n=2 threshold x0 < %.3f: x0=0.09 ->'%thr, run(2,tauv,0.09,1e7), ' x0=0.11 ->', run(2,tauv,0.11,1e7))
assert len(run(2,tauv,0.09,1e7))==1 and len(run(2,tauv,0.11,1e7))==0
