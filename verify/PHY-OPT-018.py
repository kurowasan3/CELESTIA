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
# (5) zone plate rho_j=b sqrt(j), odd zones open, plane-wave illumination: on-axis amplitude (paraxial Fresnel)
# U(q) ∝ (1/(lam q)) ∫_open exp(i pi rho^2/(lam q)) d(rho^2)  -> analytic per ring
def U_zp(qv,b,lv,K):
    u=lambda r2: np.pi*r2/(lv*qv)
    tot=0j
    for j in range(1,2*K+1,2):
        u1,u2=u((j-1)*b**2),u(j*b**2)
        tot+=(np.exp(1j*u2)-np.exp(1j*u1))/(1j)/np.pi   # = (1/(lam q)) * integral of exp(iu) d(rho^2)
    return abs(tot)
b=0.5e-3; K=20
for lv,name in [(6.5e-7,'red'),(4.0e-7,'violet')]:
    F=b**2/lv
    qs=np.linspace(F/9,F*1.6,40001); A=np.array([U_zp(qq,b,lv,K) for qq in qs])
    pk=qs[1:-1][(A[1:-1]>A[:-2])&(A[1:-1]>=A[2:])&(A[1:-1]>0.6*A.max())]
    print(name,'F=%.4f m'%F,' strong peaks at q/F =',np.round(pk/F,4))
    assert np.allclose(sorted(pk/F,reverse=True)[:4],[1,1/3,1/5,1/7],atol=2e-3)
    for m in [2,4,6]: print('   q=F/%d amplitude %.3e (vs peak %.2f)'%(m,U_zp(F/m,b,lv,K),A.max()))
print('red focus nearer than violet: F_red < F_violet ->',b**2/6.5e-7 < b**2/4.0e-7)
# exact (non-paraxial path, Rayleigh-Sommerfeld-like with 1/r) check of the main focus for a modest zone plate
from scipy.integrate import quad
def U_exact(qv,b,lv,K):
    k=2*np.pi/lv; tot=0j
    for j in range(1,2*K+1,2):
        r1,r2=b*np.sqrt(j-1),b*np.sqrt(j)
        g=lambda r: np.exp(1j*k*np.sqrt(qv**2+r**2))/np.sqrt(qv**2+r**2)*r
        tot+=quad(lambda r: g(r).real,r1,r2,limit=200)[0]+1j*quad(lambda r: g(r).imag,r1,r2,limit=200)[0]
    return abs(tot)
lv=5.0e-7; F=b**2/lv; qs=np.linspace(0.9*F,1.1*F,201); A=[U_exact(qq,b,lv,10) for qq in qs]
print('exact-path main focus at q/F =',round(qs[int(np.argmax(A))]/F,4)); assert abs(qs[int(np.argmax(A))]/F-1)<0.01
print('OK')
