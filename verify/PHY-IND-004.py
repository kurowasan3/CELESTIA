# N39-1: square coil (side a) launched away from long straight wire; escape speed
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
mu0,I,a,m,R,x,v,v0=sp.symbols('mu0 I a m R x v v0',positive=True)
kap=mu0*I*a/(2*sp.pi)
Phi=kap*sp.log((x+a)/x)
emf=sp.simplify(-sp.diff(Phi,x)*v)
i=emf/R; print('i =',sp.factor(i))
# force: sides parallel to wire, B=mu0 I/(2 pi r)
Fnet=sp.simplify(i*a*(mu0*I/(2*sp.pi*x)-mu0*I/(2*sp.pi*(x+a))))
print('F =',sp.factor(Fnet), ' equals (dPhi/dx)^2 v/R ?',sp.simplify(Fnet-sp.diff(Phi,x)**2*v/R))
J=sp.integrate(sp.diff(Phi,x)**2,(x,a,sp.oo))
vc=sp.simplify(J/(m*R)); print('vc =',vc, ' =',sp.simplify(vc/((mu0*I/(2*sp.pi))**2*a/(m*R))))
print('numeric factor 3/2-2ln2 =',float(sp.Rational(3,2)-2*sp.log(2)))
# numeric check: final speed = v0 - vc independent of v0
vals={mu0:4e-7*np.pi,I:1.0,a:1.0,m:1.0,R:1.0}
k=(2e-7)**2  # (mu0 I/2pi)^2 with I=1
# scale: use dimensionless with K=(mu0 I/2pi)^2 a^4/(mR)=1 -> choose units
K=1.0;A=1.0
def rhs(t,y):
    xx,vv=y; G=A**4/(xx**2*(xx+A)**2)
    return [vv,-K*G*vv]
vcn=K*(1.5-2*np.log(2))/A
for v0n in [3*vcn,1.5*vcn]:
    s=solve_ivp(rhs,[0,1e5],[A,v0n],rtol=1e-11,atol=1e-13)
    print('v0/vc=',v0n/vcn,' v_final/vc=',s.y[1,-1]/vcn,' (expect',v0n/vcn-1,')')
s=solve_ivp(rhs,[0,1e6],[A,0.999*vcn],rtol=1e-11,atol=1e-14)
print('just below vc: final v',s.y[1,-1],'x',s.y[0,-1])
