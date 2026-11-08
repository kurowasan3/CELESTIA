# S0403-1 rotating liquid lens: symbolic + exact ray tracing
import sympy as sp, numpy as np
from scipy.optimize import brentq
w,g,n,H,R,r,z0=sp.symbols('omega g n H R r z0',positive=True)
# (1) volume conservation
z=z0+w**2*r**2/(2*g)
V=sp.integrate(2*sp.pi*r*z,(r,0,R))
z0s=sp.solve(sp.Eq(V,sp.pi*R**2*H),z0)[0]
print('z0 =',sp.simplify(z0s))
f1=g/((n-1)*w**2); D=f1+z0s/n
print('D =',sp.simplify(D)); print('dD/d(w^2)=',sp.simplify(sp.diff(D.subs(w,sp.sqrt(sp.Symbol('W',positive=True))),sp.Symbol('W'))))
wc2=4*g*H/R**2; print('D at wc:',sp.simplify(D.subs(w,sp.sqrt(wc2))))
# numeric exact ray trace
gv,nv,Hv,Rv,wv=9.8,1.33,0.05,0.10,6.0
z0v=Hv-wv**2*Rv**2/(4*gv)
def surf(x): return z0v+wv**2*x**2/(2*gv)
def slope(x): return wv**2*x/gv
def refract(d,nrm,n1,n2):
    d=d/np.linalg.norm(d); nrm=nrm/np.linalg.norm(nrm)
    c=-np.dot(nrm,d)
    if c<0: nrm=-nrm; c=-c
    k=n1/n2; s=1-k*k*(1-c*c)
    return k*d+(k*c-np.sqrt(s))*nrm
f1v=gv/((nv-1)*wv**2)
# from below: ray up at x, exits surface; back-extend to axis
for x in [1e-3,3e-3]:
    p=np.array([x,surf(x)]); nrm=np.array([-slope(x),1.0])
    d=refract(np.array([0,1.0]),nrm,nv,1.0)
    t=-p[0]/d[0]; zc=p[1]+t*d[1]
    print('below: x=',x,' F below vertex =',z0v-zc,' theory',f1v)
# from above
for x in [1e-3,3e-3]:
    p=np.array([x,surf(x)]); nrm=np.array([-slope(x),1.0])
    d=refract(np.array([0,-1.0]),nrm,1.0,nv)
    t=-p[1]/d[1]; pb=p+t*d
    d2=refract(d,np.array([0,1.0]),nv,1.0)
    tt=-pb[0]/d2[0]; zc=pb[1]+tt*d2[1]
    print('above: x=',x,' Fprime above bottom =',zc,' theory',f1v+z0v/nv)
print('wc2',4*gv*Hv/Rv**2,' w^2',wv**2)
