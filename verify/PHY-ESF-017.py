# Q113-1: charged conducting spherical shell cut by a plane: force on cap, lift-off condition
import numpy as np, sympy as sp
th,h,R,Q,e0,sm,g=sp.symbols('theta h R Q epsilon0 sigma_m g',real=True,positive=True)
sig=Q/(4*sp.pi*R**2); p=sig**2/(2*e0)
th0=sp.acos(h/R)
Fz=sp.integrate(p*sp.cos(th)*2*sp.pi*R**2*sp.sin(th),(th,0,th0))
print('Fz',sp.factor(sp.simplify(Fz)),' vs p*pi(R^2-h^2):',sp.simplify(Fz-p*sp.pi*(R**2-h**2)))
mcap=sm*2*sp.pi*R*(R-h)
Q2=sp.solve(sp.Eq(Fz,mcap*g),Q)[0]**2; print('Q^2 threshold',sp.factor(sp.simplify(Q2)))
# limit h->R: local condition sigma^2/(2e0)=sigma_m g
print('limit check', sp.simplify((Q2/(16*sp.pi**2*R**4)/(2*e0)).subs(h,R)))
# Monte Carlo direct Coulomb force between cap (z>h) and rest, R=1, Q=1, k=1
rng=np.random.default_rng(2)
hv=0.3; N=40000
def pts(n):
    v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1)[:,None]
ests=[]
for rep in range(30):
    P=pts(N); cap=P[P[:,2]>hv]; rest=P[P[:,2]<=hv]
    q=1.0/N
    # sample subsets to limit cost
    A=cap[:3000]; B=rest[:3000]
    wA=len(cap)/len(A); wB=len(rest)/len(B)
    d=A[:,None,:]-B[None,:,:]; r=np.linalg.norm(d,axis=2)
    ests.append((q*q*wA*wB*d[:,:,2]/r**3).sum())
print('MC Fz',np.mean(ests),'+-',np.std(ests)/np.sqrt(30))
# theory with k=1/(4 pi e0)=1 -> e0=1/(4pi)
print('theory',float(Fz.subs({Q:1,R:1,h:hv,e0:1/(4*np.pi)})))
# numbers: R=0.20 m, sigma_m=0.050 kg/m^2, hemisphere h=0
vals={R:0.20,sm:0.050,g:9.8,e0:8.85e-12,h:0}
print('Q hemisphere [C]',float(sp.sqrt(Q2.subs(vals))))
