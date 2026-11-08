# Q107-1: sand on inclined belt; heat, belt work, power, inverse for mu
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,v,al,mu,lam,L=sp.symbols('g v alpha mu lambda L',positive=True)
a=g*(mu*sp.cos(al)-sp.sin(al))
t1=v/a; s1=v**2/(2*a)
f=mu*g*sp.cos(al)               # friction per unit mass
Q=f*(v*t1-s1)                   # heat per unit mass = f * relative slip
Wb=f*v*t1                       # belt work per unit mass during slip
print('t1',sp.simplify(t1),' s1',sp.simplify(s1))
print('Q',sp.simplify(Q),' Wb/Q',sp.simplify(Wb/Q))
# energy check: Wb = KE + PE(slip) + Q
print('energy check',sp.simplify(Wb-(v**2/2+g*s1*sp.sin(al)+Q)))
P=lam*(Wb+g*sp.sin(al)*(L-s1))
print('P',sp.simplify(P),' P - lam(v^2/2+gLsin+Q)=',sp.simplify(P-lam*(v**2/2+g*L*sp.sin(al)+Q)))
# inverse: Q = r v^2/2 -> mu
r=sp.symbols('r',positive=True)
sol=sp.solve(sp.Eq(sp.simplify(Q),r*v**2/2),mu); print('mu(r)=',sol)
# numeric: belt work via integration of a single grain
gv,vv,alv,muv=9.8,2.0,np.radians(20),0.8
def rhs(t,y): return [gv*(muv*np.cos(alv)-np.sin(alv))]
ev=lambda t,y:y[0]-vv; ev.terminal=True
so=solve_ivp(rhs,[0,10],[0.0],events=ev,max_step=1e-4,dense_output=True)
T=so.t_events[0][0]; ts=np.linspace(0,T,20001); u=so.sol(ts)[0]
slip=np.trapezoid(vv-u,ts); Qn=muv*gv*np.cos(alv)*slip; Wn=muv*gv*np.cos(alv)*vv*T
print('numeric Q',Qn,' formula',float(Q.subs({g:gv,v:vv,al:alv,mu:muv})),' Wb/Q',Wn/Qn)
# numbers for (4): v=2, alpha=30deg, lam=40, L=10, measured P=2240 W
vals={g:9.8,v:2.0,al:sp.pi/6,lam:40,L:10}
P0=40*(2.0**2+9.8*10*0.5); print('lam(v^2+gLsin)=',P0)
Qm=2.0**2/2+(2240-P0)/40; rr=2*Qm/2.0**2
mu4=rr*sp.tan(sp.pi/6)/(rr-1); print('Q',Qm,'mu4',float(mu4),sp.nsimplify(mu4))
print('P total check', float(P.subs(vals).subs(mu,mu4)))
print('s1 at mu4', float(s1.subs(vals).subs(mu,mu4)), ' a=',float(a.subs(vals).subs(mu,mu4)))
