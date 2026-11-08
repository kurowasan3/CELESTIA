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
# ---- (5) rewritten: neutral dipole layer still shifts the work function ----
n_,p_,tp=sp.symbols('n p tprime',positive=True)
# layer: +sigma' (outer) at z=tp, -sigma' (inner) at z=0, sigma'*tp = n p ; metal interior z<0, vacuum z>tp
sg=n_*p_/tp
def Efield(z):   # superposition of two infinite sheets (field along +z)
    return sp.Piecewise((0,z<0),(-sg/e0,z<tp),(0,True))
zz=sp.symbols('z')
dphi=-sp.integrate(-sg/e0,(zz,0,tp))      # phi(outside)-phi(inside)
print('(5a) field outside layer: 0 (sheets cancel); potential step outside-inside =',sp.simplify(dphi))
assert sp.simplify(dphi-n_*p_/e0)==0
# numerical check with a smooth (non-sheet) neutral charge distribution: rho(z) = A*(z-tp/2)*bump, dipole moment per area = n p
import numpy as np
ttn=2e-10; zg=np.linspace(-ttn,2*ttn,300001); dz=zg[1]-zg[0]
bump=np.where((zg>0)&(zg<ttn),np.sin(np.pi*zg/ttn)**2,0.0)
rho=(zg-ttn/2)*bump
dip=np.sum(rho*zg)*dz; target=1.0e18*3.0e-30
rho*=target/dip
Ez=np.cumsum(rho)*dz/eps0               # Gauss: E(z)= (1/eps0) int rho ; E(-inf)=0
phi=-np.cumsum(Ez)*dz
print('    smooth layer: net charge',np.sum(rho)*dz,' E far outside',Ez[-1],' dphi',phi[-1]-phi[0],' np/eps0',target/eps0)
assert abs(Ez[-1])<1e-6*np.abs(Ez).max() and abs((phi[-1]-phi[0])-target/eps0)<1e-3*target/eps0
WAp=WA-e*n_*p_/e0
Ep_=sp.solve(sp.Eq((WAp-WB)/e-E,0),E)[0]
print('(5b) null EMF after adsorption:',sp.simplify(Ep_))
assert sp.simplify(Ep_-((WA-WB)/e-n_*p_/e0))==0
nstar=sp.solve(Ep_,n_)[0]; print('     n for zero EMF:',nstar)
assert sp.simplify(nstar-e0*(WA-WB)/(e*p_))==0
# Kelvin current with adsorbed layer and no battery vanishes at n*
Qn=e0*S*((WAp-WB)/e)/(d0+dl*sp.sin(w*t))
print('     current at n*:',sp.simplify(sp.diff(Qn,t).subs(n_,nstar)))
# old numerical example kept as illustration: dE=0.45 V, p=3e-30 -> n
print('     example: n for 0.45 V shift =',eps0*0.45/3.0e-30,' m^-2')
print('ESF-007 (5) OK')
