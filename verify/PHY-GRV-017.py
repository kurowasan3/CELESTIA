# Q109-1: two-layer planet g(r)
import sympy as sp, numpy as np
r,R,k,G,rho=sp.symbols('r R k G rho2',positive=True)
a=R/2
Mc=sp.Rational(4,3)*sp.pi*k*rho*a**3
gm=G*(Mc+sp.Rational(4,3)*sp.pi*rho*(r**3-a**3))/r**2   # mantle a<r<R
gR=gm.subs(r,R); ga=gm.subs(r,a)
print('g(a)/g(R)=',sp.simplify(ga/gR))
print('dg/dr at R:',sp.factor(sp.simplify(sp.diff(gm,r).subs(r,R))))
rs=sp.solve(sp.diff(gm,r),r); print('stationary r',rs)
# k=3
g3=gm.subs(k,3); rstar=(4)**sp.Rational(1,3)*R/2
print('r*/R',float(rstar/R),' 1-r*/R',float(1-rstar/R))
print('gmin/gR',sp.nsimplify(sp.simplify(g3.subs(r,rstar)/g3.subs(r,R))), float(g3.subs(r,rstar)/g3.subs(r,R)))
print('g(a)/gR k=3',sp.simplify(ga.subs(k,3)/gR.subs(k,3)))
# (1) general: surface shell density rho_s, mean rho_bar: dg/dd at surface
rs_,rb,d=sp.symbols('rho_s rho_bar d',positive=True)
gd=G*(sp.Rational(4,3)*sp.pi*rb*R**3-4*sp.pi*rs_*R**2*d)/(R-d)**2
print('dg/dd at d=0:',sp.factor(sp.diff(gd,d).subs(d,0)))
# numeric scan for classification
for kv in [2,7/3,3,5,6]:
    rr=np.linspace(0.5,1,20001); gg=np.array([float(gm.subs({k:kv,R:1,G:1,rho:1,r:x})) for x in rr[::400]])
    print(kv,'monotone decreasing in r?',np.all(np.diff(gg)<=1e-12),' g(a)>g(R)?',gg[0]>gg[-1])
