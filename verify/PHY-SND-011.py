# J079-1: partial reflection standing wave of sound pressure; recover r, delta, absorption; velocity sensor
import numpy as np, sympy as sp
A=1.0; f=500.0; V=340.0; lam=V/f; k=2*np.pi/lam; r=0.5; delta=np.pi/2
w=2*np.pi*f; ls=np.linspace(0,0.6,60001); ts=np.linspace(0,1/f,400,endpoint=False)
# incident p_i(x,t)=A sin(wt-kx); reflected at x=0 lagging delta: r A sin(wt+kx-delta); point at x=-l
P=np.array([np.max(np.abs(A*np.sin(w*ts+k*l)+r*A*np.sin(w*ts-k*l-delta))) for l in ls[::50]])
lsub=ls[::50]
print('pmax/pmin',P.max()/P.min(),' expected s=',(1+r)/(1-r))
# minima positions
from scipy.signal import argrelmin
idx=argrelmin(P)[0]; print('pressure minima l:',lsub[idx], ' spacing',np.diff(lsub[idx]),' lam/2=',lam/2)
# formula: delta = pi(1-2 l0/D) mod 2pi
D=lam/2; l0=lsub[idx][0]; print('l0=',l0,' delta from formula',np.pi*(1-2*l0/D),' true',delta)
# given data s=3, D=0.34, l0=0.085, f=500
s=3; print('V=',2*0.34*500,' r=',(s-1)/(s+1),' absorbed',4*s/(s+1)**2,' delta/pi',(1-2*0.085/0.34))
# velocity: u ∝ p_i - p_r
U=np.array([np.max(np.abs(A*np.sin(w*ts+k*l)-r*A*np.sin(w*ts-k*l-delta))) for l in lsub])
iu=argrelmin(U)[0]; print('velocity minima',lsub[iu],' ratio',U.max()/U.min())
# symbolic amplitude
l,kk,rr,dd=sp.symbols('l k r delta',real=True)
z=sp.exp(sp.I*kk*l)+rr*sp.exp(sp.I*(-kk*l-dd))
print(sp.simplify(sp.expand(z*sp.conjugate(z)).rewrite(sp.cos)))
