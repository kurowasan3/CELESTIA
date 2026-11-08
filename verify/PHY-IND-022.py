# T20152-1: two rods (m, 2m) on inclined rails in vertical B
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,g,R,B,L,th,t,u=sp.symbols('m g R B L theta t u',positive=True)
c=sp.cos(th); s=sp.sin(th)
# rod1 alone (rod2 fixed): current = B L c v /(2R), force along rail B I L c
U=sp.solve(sp.Eq(m*g*s, B**2*L**2*c**2*u/(2*R)),u)[0]; print('(1) u=',U)
v1,v2=sp.Function('v1')(t),sp.Function('v2')(t)
I=B*L*c*(v1-v2)/(2*R)
eqs=[sp.Eq(m*v1.diff(t), m*g*s - B*I*L*c), sp.Eq(2*m*v2.diff(t), 2*m*g*s + B*I*L*c)]
sol=sp.dsolve(eqs,[v1,v2],ics={v1.subs(t,0):U, v2.subs(t,0):0})
V1=sp.simplify(sol[0].rhs); V2=sp.simplify(sol[1].rhs)
tau=4*m*R/(3*B**2*L**2*c**2)
V1c=U/3+g*s*t+2*U/3*sp.exp(-t/tau); V2c=U/3+g*s*t-U/3*sp.exp(-t/tau)
print('v1 check',sp.simplify(V1-V1c),' v2 check',sp.simplify(V2-V2c))
print('(2) a1(0)',sp.simplify(V1.diff(t).subs(t,0)),' a2(0)',sp.simplify(V2.diff(t).subs(t,0)))
ds=sp.integrate(sp.simplify(V1-V2),(t,0,sp.oo)); print('(5) ds=',sp.simplify(ds), ' expected', sp.simplify(U*tau))
Q=sp.integrate(sp.simplify(I.subs({v1:V1,v2:V2})**2*2*R),(t,0,sp.oo)); print('(6) Q=',sp.simplify(Q),' expected m u^2/3 =',sp.simplify(m*U**2/3))
# energy check: total energy balance
# numeric integration incl. phase before release
vals={m:0.1,g:9.8,R:0.5,B:0.8,L:0.5,th:0.4}
mm,gg,RR,BB,LL,tt=[vals[k] for k in (m,g,R,B,L,th)]
k=BB**2*LL**2*np.cos(tt)**2/(2*RR)
def f(t,y):
    v1,v2,x1,x2=y; F=k*(v1-v2)
    return [gg*np.sin(tt)-F/mm, gg*np.sin(tt)+F/(2*mm), v1, v2]
uu=float(U.subs(vals)); tauv=float(tau.subs(vals))
s=solve_ivp(f,[0,40*tauv],[uu,0,0,0],rtol=1e-11,atol=1e-12)
print('num ds',s.y[2,-1]-s.y[3,-1],' theory',float((U*tau).subs(vals)))
