# K20002-1 verification: surface dipole layer, contact potential, Kelvin probe
import sympy as sp, numpy as np
e0,S,d0,dl,w,t,WA,WB,e,E=sp.symbols('epsilon0 S d0 delta omega t W_A W_B e E',positive=True)
# (1) numeric sigma
eps0=8.85e-12; ech=1.60e-19
W=4.5*ech; tt=1.0e-10
sig=eps0*W/(ech*tt); print('(1) sigma=',sig,'C/m^2  electrons per m^2=',sig/ech)
# (2)(3) charge with battery E (positive terminal on A side): gap voltage = (WA-WB)/e - E
Vg=(WA-WB)/e-E
Q=e0*S*Vg/(d0+dl*sp.sin(w*t))
I=sp.diff(Q,t)
I1=sp.series(I.subs(dl,sp.Symbol('k')*d0),sp.Symbol('k'),0,2).removeO()
print('(3) I to first order in delta:',sp.simplify(I1.subs(sp.Symbol('k'),dl/d0)))
print('    null E =',sp.solve(Vg,E))
# numerical check of amplitude by finite sampling
vals={e0:eps0,S:1e-4,d0:1e-3,dl:1e-5,w:2*np.pi*500,WA:5.1*ech,WB:4.3*ech,e:ech,E:0}
f=sp.lambdify(t,I.subs(vals)); ts=np.linspace(0,2e-3,20001)
print('    num amp',np.max(np.abs(f(ts))),' formula',float((e0*S*(WA-WB)/e*dl*w/d0**2).subs(vals)))
# (5) adsorption
p=3.0e-30; dE=0.62-0.17
n=eps0*dE/p; print('(5) n=',n)
print('    n* for null without battery (WA-WB=0.62 eV):',eps0*0.62/p)
