# N59-1 X-ray steps + electron Bragg with inner potential
import numpy as np, sympy as sp
h=6.6e-34; e=1.6e-19; c=3.0e8; m=9.0e-31
th=np.radians(30); s=np.sin(th)
# (1) thresholds V_n = n*h*c/(2 e d sin) ; given step 6.2e3 V
dV=6.2e3; d=h*c/(2*e*s*dV); print('(1) d =',d)
d=2.0e-10
# (2) at 20 kV: lam_min, reflected wavelengths d*2s/n >= lam_min
lmin=h*c/(e*20e3); print('(2) lam_min',lmin, [2*d*s/n for n in range(1,10) if 2*d*s/n>=lmin])
# (3) symbolic electron condition with refraction (independent: wave-vector approach)
V,V1,n,mm,ee,hh,dd,t=sp.symbols('V V1 n m e h d theta',positive=True)
p=sp.sqrt(2*mm*ee*V); pin=sp.sqrt(2*mm*ee*(V+V1))
pz_in=sp.sqrt(pin**2-(p*sp.cos(t))**2)    # parallel component conserved
cond=sp.Eq(2*dd*pz_in, n*hh)                 # path diff 2 d sin(phi) = n lambda'
Vsol=sp.solve(cond,V)[0]; print('(3) V_n =',sp.simplify(Vsol))
K=h**2/(8*m*e*d**2); print('K =',K)
Vn=lambda n,V1: (K*n**2-V1)/s**2
print('peaks V1=15:',[(k,round(Vn(k,15),1)) for k in range(4,10)])
# brute-force numeric: wave in crystal, compute phase difference for scanned V
def phase_ok(V,V1):
    lam_in=h/np.sqrt(2*m*e*(V+V1)); pin=np.sqrt(2*m*e*(V+V1)); p=np.sqrt(2*m*e*V)
    phi=np.arccos(p*np.cos(th)/pin); return 2*d*np.sin(phi)/lam_in
Vs=np.linspace(1000,2000,200001); r=phase_ok(Vs,15.3)
print('integer crossings:',Vs[np.where(np.diff(np.floor(r))!=0)[0]])

# ---- new (4): shift -V1/sin^2(theta), direction, missing low orders ----
Vsol0=Vsol.subs(V1,0)
print('V_n(V1) - V_n(0) =',sp.simplify(Vsol-Vsol0),'  expected -V1/sin^2:',sp.simplify(Vsol-Vsol0+V1/sp.sin(t)**2))
# numeric brute force: peak positions with and without V1 at two angles
for thdeg in [30,10]:
    thr=np.radians(thdeg)
    def phase(V,V1):
        p=np.sqrt(2*m*e*V); pin=np.sqrt(2*m*e*(V+V1)); return 2*d*np.sqrt(pin**2-(p*np.cos(thr))**2)/h
    Vs2=np.linspace(1,20000,2000001)
    def crossings(V1): r=phase(Vs2,V1); return Vs2[np.where(np.diff(np.floor(r))!=0)[0]][:3]
    c0=crossings(0.0); c15=crossings(15.0)   # with V1=15 order n=1 is absent, so c15[0] is n=2
    print(' theta=%d: shifts for n=2,3:'%thdeg,c15[0]-c0[1],c15[1]-c0[2],' predicted',-15/np.sin(thr)**2)
# missing orders: n^2 K <= V1 -> no V>0 satisfies
print('K=',K,' orders missing for V1=15:',[k for k in range(1,5) if k**2*K<=15])
print('min internal normal momentum*2d/h at V->0 (=max missing n):',2*d*np.sqrt(2*m*e*15)/h)
