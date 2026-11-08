# N57-1 Compton: successive scatterings
import sympy as sp
E0,E,m,c,th,a,N=sp.symbols('E0 E m c theta alpha N',positive=True)
mc2=sp.Symbol('M',positive=True)  # mc^2
def scat(Ein,ang): return 1/(1/Ein+(1-sp.cos(ang))/mc2)
# (1) max KE
K=sp.simplify(E0-scat(E0,sp.pi)); print('(1) Kmax =',sp.factor(K))
# (2) 90deg photon = E0/3
sol=sp.solve(sp.Eq(scat(E0,sp.pi/2),E0/3),E0); print('(2) E0 =',sol)
E0v=2*mc2
# (3) recoil angle via momentum conservation, independent: solve vector eq
Ep=scat(E0v,sp.pi/2); px=E0v; py=Ep  # electron momentum*c components (photon goes -y)
print('(3) tan phi =',sp.simplify(py/px))
# full relativistic check of Compton formula (independent): energy & momentum conservation
Ein,Eout,phi=sp.symbols('Ein Eout phi',positive=True)
pe2=(Ein**2+Eout**2-2*Ein*Eout*sp.cos(th))  # (pc)^2
poly=sp.expand((Ein+mc2-Eout)**2-(mc2**2+pe2))
s=sp.solve(poly,Eout); print('rel check 1/Eout-1/Ein =',[sp.simplify(1/x-1/Ein) for x in s])
# (4) two steps alpha vs single 2alpha
d=sp.simplify(2*(1-sp.cos(a))-(1-sp.cos(2*a))); print('(4) diff 1/E (two-step - single)*mc2 =',sp.factor(sp.expand_trig(d)))
# (5) N steps pi/N, E0=2mc2
EN=1/(1/E0v+N*(1-sp.cos(sp.pi/N))/mc2)
for n in [1,2,3,4,6,10,100]: print(n, float(EN.subs(N,n)/mc2))
print('limit', sp.limit(EN/mc2,N,sp.oo))
# numeric sequential simulation with explicit momentum bookkeeping (independent): track photon through N electrons
import numpy as np
def sim(n):
    Eph=2.0
    for k in range(n):
        al=np.pi/n
        # solve energy-momentum exact relativistic numerically
        from scipy.optimize import brentq
        f=lambda Eo: (Eph+1-Eo)-np.sqrt(1+Eph**2+Eo**2-2*Eph*Eo*np.cos(al))
        Eph=brentq(f,1e-9,Eph)
    return Eph
print([round(sim(n),6) for n in [1,2,3,4,6,10,100]])
