# M64-1 verification: string resonance with buoyancy-reduced tension
import sympy as sp
from fractions import Fraction as F
a1,a2=0.60,0.75
sig,h,s,rho,g=2.7e3,0.10,2.0e-4,1.0e-3,9.8
lam=2*(a2-a1); n2=round(a2/(lam/2)); n1=round(a1/(lam/2))
print("lam",lam,"n1,n2",n1,n2)
T0=sig*s*h*g; v=(T0/rho)**0.5; f=v/lam
print("T0",T0,"v",v,"f",f)
# (2) exactly one resonance (n=6) while immersing 0<d<=h: T(d)=T0(1-rhoL d/(sig h)); resonance when T=T0(5/n)^2
# need (5/6)^2 >= 1-rhoL/sig > (5/7)^2
lo=sig*(1-sp.Rational(25,36)); hi=sig*(1-sp.Rational(25,49))
print("rhoL range",float(lo),float(hi), sp.nsimplify(1-sp.Rational(25,36)), 1-sp.Rational(25,49))
# brute check by scanning rhoL
import numpy as np
def count(rhoL):
    Tmin=T0*(1-rhoL/sig); c=0
    for n in range(6,100):
        if T0*(5/n)**2>=Tmin-1e-12: c+=1
    return c
for r in [800,825,900,1000,1300,1322,1323,1400]:
    print(r,count(r))
# (3) water
rho0=1.0e3
d=sig*h*(1-25/36)/rho0
print("d",d)
# (4) full immersion
lam2=lam*(1-rho0/sig)**0.5
ks=[k*lam2/2 for k in range(1,12)]
best=min(ks,key=lambda L:abs(L-a2)); print("lam2",lam2,"best",best,"shift",best-a2)
print("right option",[L for L in ks if L>a2][0]-a2)
# independent: direct frequency matching f = k/(2a) sqrt(T/rho)
T=T0*(1-rho0/sig)
for k in range(5,9):
    a=k/(2*f)*(T/rho)**0.5; print(k,a)
