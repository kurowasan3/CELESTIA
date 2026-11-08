# M45-1: diametral tunnel through uniform earth; released from rest at r=2R above entrance.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,R,r=sp.symbols('g R r',positive=True)
GM=g*R**2
vs=sp.sqrt(2*(GM/R-GM/(2*R))); print('(1) v at surface',sp.simplify(vs))
w=sp.sqrt(g/R)
A=sp.sqrt(R**2+vs**2/w**2); print('(2) amplitude',sp.simplify(A),' vmax',sp.simplify(A*w))
tin=2*sp.asin(R/A)/w; print('(3) t_in',sp.simplify(tin))
# outside: degenerate ellipse a=R (from 2R to focus O), period 2pi sqrt(a^3/GM)
Tell=2*sp.pi*sp.sqrt(R**3/GM)
frac=(sp.pi*R**2/4+R**2/2)/(sp.pi*R**2)
tout=sp.simplify(Tell*frac); print('(4) t_out',tout)
# direct integral check: t = int_R^{2R} dr / sqrt(2GM(1/r-1/(2R)))
tint=sp.integrate(1/sp.sqrt(2*GM*(1/r-1/(2*R))),(r,R,2*R))
print('    integral',sp.simplify(tint), ' diff',sp.N((tint-tout).subs({g:9.8,R:1})))
P=2*(2*tout+tin); print('(5) period',sp.simplify(P))
# numeric
gg,RR=9.8,6.4e6
def f(t,y):
    x=y[0]; a=-gg*x/RR if abs(x)<=RR else -gg*RR**2*np.sign(x)/x**2
    return [y[1],a]
ev=lambda t,y: y[1]; ev.direction=1
r_=solve_ivp(f,[0,30000],[2*RR,0],rtol=1e-12,atol=1e-6,max_step=2,events=ev)
print('num period',r_.t_events[0][:2],' th',float(P.subs({g:gg,R:RR})))
print('num vmax',np.abs(r_.y[1]).max(),' th',float((A*w).subs({g:gg,R:RR})))
