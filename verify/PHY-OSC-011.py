# M39-1: Q(2m) held, P(m) hangs by rubber band (k, natural L). Release Q.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,k,g,L,t=sp.symbols('m k g L t',positive=True)
x0=m*g/k
mu=2*m*m/(3*m); w=sp.sqrt(k/mu)
print('omega^2',sp.simplify(w**2))
# half-band from Q to CM: natural length L/3 -> constant 3k ; Q mass 2m
print('half-band check', sp.simplify(3*k/(2*m)-w**2))
aQ=(2*m*g+k*x0)/(2*m); aP=(m*g-k*x0)/m
print('(1) aQ,aP',sp.simplify(aQ),sp.simplify(aP))
t1=sp.pi/2/w; u=x0*w
print('(2) t1',sp.simplify(t1),' (3) u',sp.simplify(u))
vcm=g*t1
print('(3) vQ',sp.simplify(vcm+u/3),' vP',sp.simplify(vcm-2*u/3))
T=t1+L/u
dQ=g*T**2/2+(L+x0)/3
print('(4) T',sp.simplify(T),' dQ',sp.expand(sp.simplify(dQ)))
vP=g*t - sp.Rational(2,3)*x0*w*sp.sin(w*t)
print('(5) vP(t)=',sp.simplify(vP), ' = g(t-sin(wt)/w)?', sp.simplify(vP-g*(t-sp.sin(w*t)/w)))
print('    vP at T >0 ; g t1-2u/3 =',sp.simplify(g*t1-2*u/3))
# numeric
mm,kk,gg,LL=0.2,40.0,9.8,0.5
def f(tt,y):
    yQ,vQ,yP,vP_=y   # downward coordinates
    s=(yP-yQ)-LL
    T_=kk*s if s>0 else 0.0
    return [vQ,gg+T_/(2*mm),vP_,gg-T_/mm]
ev=lambda tt,y: y[2]-y[0]; ev.terminal=True; ev.direction=-1
x0n=mm*gg/kk
r=solve_ivp(f,[0,5],[0,0,LL+x0n,0],events=ev,rtol=1e-11,atol=1e-13,max_step=1e-4,dense_output=True)
Tn=r.t_events[0][0]; yQ,vQ,yP,vPn=r.y_events[0][0]
sub={m:mm,k:kk,g:gg,L:LL}
print('num T',Tn,' th',float(T.subs(sub)))
print('num Q fall',yQ,' th',float(dQ.subs(sub)))
print('num vQ,vP',vQ,vPn,' th', float((g*T+u/3).subs(sub)), float((g*T-2*u/3).subs(sub)))
ts=np.linspace(0,Tn,4000); vp=r.sol(ts)[3]
print('min vP over motion', vp.min())
