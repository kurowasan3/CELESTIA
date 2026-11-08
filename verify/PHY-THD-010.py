import sympy as sp, numpy as np
from scipy.integrate import quad
R,T0,p0,V0,S=sp.symbols('R T0 p0 V0 S',positive=True)
g=sp.Rational(5,3)
VB=V0-2*(V0/4); TB=T0*(V0/VB)**(g-1); pB=p0*(V0/VB)**g
pA=2*pB-p0; TA=sp.simplify(pA*sp.Rational(5,4)*V0/(p0*V0/T0))   # R = p0 V0/T0
print('VB',VB,'TB',sp.nsimplify(TB/T0),sp.N(TB/T0),'pB/p0',sp.N(pB/p0),'pA/p0',sp.N(pA/p0),'TA/T0',sp.simplify(TA/T0),sp.N(TA/T0))
W0=p0*V0/4
q=sp.Rational(3,2)*(p0*V0/T0)*(TA-T0)+sp.Rational(3,2)*(p0*V0/T0)*(TB-T0)-W0
print('q/(RT0)=',sp.nsimplify(sp.simplify(q/(p0*V0))),sp.N(q/(p0*V0)))
# numeric independent: integrate work done by A along quasi-static path, units p0=V0=T0=R=1
def pA_of(u):  # u = S l / V0
    vb=1-2*u; pb=vb**(-5/3); return 2*pb-1
WA,_=quad(pA_of,0,0.25)
dUA=1.5*(pA_of(0.25)*1.25-1)
print('numeric WA',WA,'q',dUA+WA)
assert abs(dUA+WA-float(q/(p0*V0)))<1e-9
# sign checks: pA>0 along path, and heat into A positive along path
us=np.linspace(0,0.25,2001); 
dq=[1.5*np.gradient([pA_of(x)*(1+x) for x in us],us)[i]+pA_of(us[i]) for i in range(len(us))]
print('min dq/du',min(dq)); print('OK')
