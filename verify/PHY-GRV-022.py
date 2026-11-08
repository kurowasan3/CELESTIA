# S0611-1: light in an accelerating box, engine cut-off, clock lag; numeric crest tracking in inertial frame S
import numpy as np, sympy as sp
from scipy.optimize import brentq
c=1.0; h=1.0; a=2e-4; T=30.0; f0=200.0
def v(t): return 0.0 if t<0 else (a*t if t<T else a*T)
def y(t):  # floor position
    if t<0: return 0.0
    if t<T: return 0.5*a*t*t
    return 0.5*a*T*T + a*T*(t-T)
def t_recv(te):
    ye=y(te)
    return brentq(lambda tr: ye + c*(tr-te) - (y(tr)+h), te, te+10*h/c)
# crests emitted at A at times k/f0 (source frequency f0 in S; nonrelativistic)
tes=np.arange(-5.0, T+5.0, 1/f0)
trs=np.array([t_recv(te) for te in tes])
fB=(1/f0)/np.diff(trs)*f0  # received frequency
tm=0.5*(trs[1:]+trs[:-1])
def model(tr):
    if tr<0: return 1
    if tr<h/c: return 1-a*tr/c
    if tr<T: return 1-a*h/c**2
    if tr<T+h/c: return 1-a*(T+h/c-tr)/c
    return 1
err=max(abs(fB[i]/f0-model(tm[i])) for i in range(len(tm)))
print("max |fB/f0 - trapezoid model| =",err," (shift size a h/c^2 =",a*h/c**2,")")
# plateau value
mid=np.argmin(abs(tm-T/2)); print("plateau fB/f0-1 =",fB[mid]/f0-1," model -ah/c^2 =",-a*h/c**2)
# clock lag: lag = tr - te (B reads A's clock te at its time tr)
print("lag before:",trs[0]-tes[0]," lag after:",trs[-1]-tes[-1]," increase:",(trs[-1]-tes[-1])-(trs[0]-tes[0])," model a h T/c^2 =",a*h*T/c**2)
# exact lag after: h/(c-aT)
print("exact h/(c-aT)-h/c =",h/(c-a*T)-h/c)
# symbolic first-order check of free-fall cancellation in Earth frame
g,hh,cc=sp.symbols('g h c',positive=True)
expr=(1-g*hh/cc**2)*(1+(g*hh/cc)/cc)
print(sp.series(expr.subs(g,sp.Symbol('e')*g),sp.Symbol('e'),0,2))
