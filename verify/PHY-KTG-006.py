import numpy as np, sympy as sp
kB=1.38e-23; c=3.0e8
# Monte Carlo: Maxwell vx, observed lambda = lam0(1 - vx/c); rms width vs formula
rng=np.random.default_rng(1)
m=3.8e-26; T=500; lam0=5.0e-7
vx=rng.normal(0,np.sqrt(kB*T/m),2_000_000)
lam=lam0*(1-vx/c)
print('MC rms dlam',np.std(lam),' formula',lam0*np.sqrt(kB*T/(m*c**2)))
r=np.sqrt(kB*T/(m*c**2)); print('r =',r)
# grating: z_k = L tan(asin(k lam/d)); spread ratio small-angle
L=1.0; d=1.0e-5
for k in (1,2,3):
    z=L*np.tan(np.arcsin(k*lam/d)); print(k,'z mean',z.mean(),'sigma/z',z.std()/z.mean())
print('sigma_2 =',2*L*lam0*r/d)
# adiabatic V->V/8 monatomic: T*8^(2/3)=4, r*2
print('ratio after compress', np.sqrt(8**(2/3)))
# free expansion to 2V (T same) then quasi-static back to V: T*2^(2/3) -> r*2^(1/3)
print('free exp then compress', 2**(1/3))
# T from r: T = m c^2 r^2/kB
print('T back', m*c**2*r**2/kB)
