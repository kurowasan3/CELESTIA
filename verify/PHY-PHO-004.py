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
# (4) inverse: given 1.30e3, 1.79e3 consecutive, only two in 1.0-2.0 kV
Va,Vb=1.30e3,1.79e3
tn=(Vb-Va)*s**2/K; print('2n+1 =',tn)
nn=round((tn-1)/2); V1v=K*nn**2-Va*s**2; print('n=',nn,'V1=',V1v, 'check via Vb:',K*(nn+1)**2-Vb*s**2)
print('in-range peaks with V1=%.1f:'%V1v,[round(Vn(k,V1v)) for k in range(1,12) if 1000<=Vn(k,V1v)<=2000])
# alternative n candidates consistency
for k in range(4,9):
    v1a=K*k**2-Va*s**2; v1b=K*(k+1)**2-Vb*s**2; print(k,round(v1a,1),round(v1b,1))
# brute-force numeric: wave in crystal, compute phase difference for scanned V
def phase_ok(V,V1):
    lam_in=h/np.sqrt(2*m*e*(V+V1)); pin=np.sqrt(2*m*e*(V+V1)); p=np.sqrt(2*m*e*V)
    phi=np.arccos(p*np.cos(th)/pin); return 2*d*np.sin(phi)/lam_in
Vs=np.linspace(1000,2000,200001); r=phase_ok(Vs,15.3)
print('integer crossings:',Vs[np.where(np.diff(np.floor(r))!=0)[0]])
