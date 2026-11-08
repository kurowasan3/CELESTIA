# P48-1 verification: source moving toward wall, reverses at x=b; observer at x=0, wall at x=L
import numpy as np, sympy as sp
from scipy.optimize import brentq
f,c,v=sp.symbols('f c v',positive=True)
fd=f*c/(c+v); fr=f*c/(c-v)
B=sp.simplify(fr-fd); print('beat before',B)  # 2cfv/(c^2-v^2)
B2=sp.simplify(f*c/(c-v)-f*c/(c+v)); print('beat after',B2, sp.simplify(B-B2))
Bstop=sp.simplify(fr-f); print('stop case beat', Bstop)
# numeric retarded-time simulation
fv,cv,vv,b,L=680.0,340.0,1.0,34.0,119.0
tr=b/vv
def X(t): return vv*t if t<=tr else b-vv*(t-tr)
def phase_direct(t):   # emission time te: t = te + X(te)/c
    te=brentq(lambda s: s+X(s)/cv-t,0,t); return 2*np.pi*fv*te
def phase_refl(t):     # path X->wall->0 : (L-X)+L
    te=brentq(lambda s: s+(2*L-X(s))/cv-t,0,t); return 2*np.pi*fv*te
def freq(ph,t,h=1e-3): return (ph(t+h)-ph(t-h))/(2*h)/(2*np.pi)
for t in [20,34.05,34.2,34.4,34.55,34.7,40]:
    fd_,fr_=freq(phase_direct,t),freq(phase_refl,t)
    print(t, round(fd_,4), round(fr_,4), 'beat',round(abs(fr_-fd_),4))
T1=b/vv+b/cv; T2=b/vv+(2*L-b)/cv; print('T1,T2',T1,T2,'tau',T2-T1)
# inverse
Bm,T1m,T2m=4.0,34.1,34.6
vs=sp.symbols('vs',positive=True)
vsol=sp.nsolve(2*fv*cv*vs/(cv**2-vs**2)-Bm,vs,1); print('v',vsol)
bsol=vsol*cv*T1m/(cv+vsol); Lsol=(cv*T2m-bsol*(cv-vsol)/vsol)/2; print('b',bsol,'L',Lsol)
# L-b from tau
print('L-b', cv*(T2m-T1m)/2)
