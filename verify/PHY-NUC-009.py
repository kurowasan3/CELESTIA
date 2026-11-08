# N69-1: recoil-distance method with unknown zero offset; Doppler at 60 deg
import numpy as np, sympy as sp
c=3.0e8; E0=1.4340; E60=1.4388
v=c*(E60/E0-1)/np.cos(np.pi/3); print('v=',v)
# exact formula E = E0/(1 - v cos/c)
print('v exact-formula',c*(1-E0/E60)/0.5)
x1,x2,f1,f2=0.10e-3,0.26e-3,0.25,0.0625
x0,L=sp.symbols('x0 L',positive=True)   # L = vT
s=sp.solve([sp.Eq(sp.Rational(1,4),2**(-(sp.Rational(1,10)+x0)/L)), sp.Eq(sp.Rational(1,16),2**(-(sp.Rational(26,100)+x0)/L))],[x0,L],dict=True)
print(s)
L=0.08e-3; print('T=',L/2.0e6, ' with v computed',L/v)
print('x0 = 0.06 mm; reading for 50%:',(L-0.06e-3)*1e3,'mm')
# Monte Carlo
rng=np.random.default_rng(1); T=4.0e-11; vv=2.0e6; tau=T/np.log(2)
t=rng.exponential(tau,1000000)
for x in [0.10e-3,0.26e-3,0.02e-3]:
    d=x+0.06e-3; print(x*1e3,'mm -> unshifted frac',np.mean(t>d/vv))
