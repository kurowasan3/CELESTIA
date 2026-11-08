# S0310-1: balloon drifting with wind; ground observer. Ground-frame brute force.
import numpy as np, sympy as sp
from scipy.optimize import brentq
V=340.0; w=20.0; D=500.0; f0=500.0
O=np.array([0.0,D])
def arrival(te):
    P=np.array([w*te,0.0])          # balloon position at emission (passes H=(0,0) at t=0)
    # in ground frame the wavefront circle centre drifts with wind: centre P+w*tau, radius V*tau
    g=lambda tau: np.linalg.norm(O-P-np.array([w*tau,0]))-V*tau
    tau=brentq(g,1e-9,1e4)
    return te+tau
def freq(te,h=1e-4):
    return f0*2*h/(arrival(te+h)-arrival(te-h))
# (1) far upwind / far downwind
print('(1) far approach', freq(-2e3), ' theory', f0*(V+w)/V)
print('    far recede ', freq(2e3), ' theory', f0*(V-w)/V)
# (2) emitted at H
print('(2) at H', freq(0.0), ' theory', f0*(1-w**2/V**2))
# (3) when heard frequency equals f0
te0=brentq(lambda te: freq(te)-f0,-200,0)
print('(3) emission x =',w*te0,' theory', -D*w/V, ' balloon position at hearing =', w*arrival(te0),' (0 expected)')
# general rule: f = f0(1 + w.n/V), n from balloon's position at reception to O
for te in [-30,-5,3,17]:
    ta=arrival(te); B=np.array([w*ta,0]); n=(O-B)/np.linalg.norm(O-B)
    print('    rule check te=%g: %.6f vs %.6f'%(te,freq(te),f0*(1+w*n[0]/V)))
# (4) harmonic mean
fp=f0*(V+w)/V; fm=f0*(V-w)/V
print('(4)', 2*fp*fm/(fp+fm), f0*(1-w**2/V**2))
a,b=sp.symbols('f_p f_m',positive=True)
F0=(a+b)/2; r=(a-b)/(a+b)
print('    symbolic', sp.simplify(F0*(1-r**2) - 2*a*b/(a+b)))
