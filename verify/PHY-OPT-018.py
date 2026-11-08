# P55-1: circular hole, point source P at p, axis point Q at q; Fresnel zones
import sympy as sp, numpy as np
from scipy.integrate import quad
rho,p,q,lam,a=sp.symbols('rho p q lambda a',positive=True)
Lpath=sp.sqrt(p**2+rho**2)+sp.sqrt(q**2+rho**2)-(p+q)
print('path excess ~',sp.simplify(sp.series(Lpath,rho,0,4).removeO()))
# zone boundary rho_j^2 = j*lam*pq/(p+q) -> area per zone
print('zone area =',sp.simplify(sp.pi*lam*p*q/(p+q)))
N=a**2*(1/p+1/q)/lam
av,lv=1.0e-3,5.0e-7
f=av**2/lv; print('a^2/lam =',f)
for pv in [4.0,1.6]:
    for Nt in [1,2,3]:
        inv=Nt*lv/av**2-1/pv
        print(' p=',pv,' N=',Nt,' q=', (1/inv if inv>0 else 'none'))
    print(' N(q->inf)=',av**2/(lv*pv))
# exact on-axis Fresnel-Kirchhoff (scalar, obliquity~1) amplitude: integral over hole
def amp(pv,qv):
    k=2*np.pi/lv
    g=lambda r: np.exp(1j*k*(np.sqrt(pv**2+r**2)+np.sqrt(qv**2+r**2)))/(np.sqrt(pv**2+r**2)*np.sqrt(qv**2+r**2))*r
    re=quad(lambda r: g(r).real,0,av,limit=2000)[0]; im=quad(lambda r: g(r).imag,0,av,limit=2000)[0]
    return abs(re+1j*im)**2*(pv+qv)**2  # normalised-ish
for pv in [4.0,1.6]:
    qs=np.linspace(0.6,12,6000); I=np.array([amp(pv,qq) for qq in qs])
    # scan from far to near: first local extremum
    rev=I[::-1]; qr=qs[::-1]
    ext=[]
    for i in range(1,len(rev)-1):
        if rev[i]>rev[i-1] and rev[i]>=rev[i+1]: ext.append(('max',round(qr[i],3)))
        if rev[i]<rev[i-1] and rev[i]<=rev[i+1]: ext.append(('min',round(qr[i],3),round(rev[i]/I.max(),4)))
    print('p=',pv,'extrema from far:',ext[:4])
# model: I ~ sin^2(pi N /2)
print('model I(N)=4 sin^2(pi N/2) relative to free propagation; N=1->4, N=2->0')
