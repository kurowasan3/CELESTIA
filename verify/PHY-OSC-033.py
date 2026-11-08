# J035-1: light ball B (m) hits A (M) on spring; second collision condition
import sympy as sp, numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
m,M,v0,k=sp.symbols('m M v0 k',positive=True)
vA,vB=sp.symbols('vA vB')
s=sp.solve([sp.Eq(m*(-v0),m*vB+M*vA), sp.Eq(vB-vA, v0)],[vA,vB])  # B moving left (-v0); elastic: relative velocity reversed
print(s)   # vA negative(left), vB positive(right) for M>m
r=sp.simplify(s[vB]/(-s[vA])); print('r=',r)
x0=brentq(lambda x: np.tan(x)-x, np.pi+1e-6, 1.5*np.pi-1e-6); rc=-np.cos(x0)
print('x0',x0,'rc',rc,'M/m crit',1+2*rc, 'sin x0', np.sin(x0))
# numerical check: simulate positions for several mass ratios; find whether gap closes
def second(Mm,v0=1.0,kk=1.0,m=1.0):
    M=Mm*m; w=np.sqrt(kk/M); vA=-2*m*v0/(m+M); vB=(M-m)*v0/(m+M)
    t=np.linspace(1e-6,40/w,400001); xA=vA/w*np.sin(w*t); xB=vB*t
    gap=xB-xA; i=np.argmax(gap<=0) if np.any(gap<=0) else None
    return None if i is None else (w*t[i], xA[i])
for Mm in [1.0,1.2,1.4,1.43,1.44,1.5,2]:
    print(Mm, second(Mm))
# M=m exact: phase pi, at x=0
print('M=m check:',second(1.0))
