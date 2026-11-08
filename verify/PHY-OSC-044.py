import sympy as sp, numpy as np
a,h,g,w,f,r,r0,z0=sp.symbols('a h g omega f r r0 z0',positive=True)
k=w**2/(2*g)
# regime 1: center above bottom: z=z0+k r^2
V1=sp.integrate(2*sp.pi*r*(z0+k*r**2),(r,0,a))
z0s=sp.solve(sp.Eq(V1,sp.pi*a**2*f*h),z0)[0]
W1=sp.solve(sp.Eq(z0s+k*a**2,h),w)
print('regime1 Omega',W1,' center height at Omega',sp.simplify(z0s.subs(w,W1[0])))
# regime 2: dry spot radius r0: z=k(r^2-r0^2)
V2=sp.integrate(2*sp.pi*r*k*(r**2-r0**2),(r,r0,a))
sol=sp.solve([sp.Eq(V2,sp.pi*a**2*f*h),sp.Eq(k*(a**2-r0**2),h)],[w,r0],dict=True)
print('regime2',sol)
# numeric brute: for f=1/4 and f=2/3, find max omega with volume conserved
def rim(om,fv,av=1,hv=1,gv=9.8):
    kk=om**2/(2*gv)
    # try regime1
    z0v=fv*hv-kk*av**2/2
    if z0v>=0: return z0v+kk*av**2
    # regime2: solve pi*kk*(a^2-r0^2)^2/2 = pi a^2 f h
    X=np.sqrt(2*av**2*fv*hv/kk)  # a^2-r0^2
    return kk*X
from scipy.optimize import brentq
for fv in [0.25,0.4,0.5,2/3]:
    Om=brentq(lambda o:rim(o,fv)-1,1e-3,100)
    pred= 2*np.sqrt(9.8*(1-fv)) if fv>=0.5 else np.sqrt(9.8/fv)
    print(fv,Om,pred)
