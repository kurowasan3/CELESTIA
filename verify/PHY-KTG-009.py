import sympy as sp, numpy as np
r,x,u,m,k=sp.symbols('r x u m k',positive=True)
# container: N1 (speed u), N2=x N1 (speed 2u). flux ∝ n v  -> R2/R1 = 2x
xs=sp.solve(sp.Eq(r,2*x),x)[0]; print('x=',xs)
v2gas=(u**2+x*4*u**2)/(1+x)
T=sp.simplify(m*v2gas/(3*k)).subs(x,xs); print('T=',sp.factor(T))
v2beam=(u**2+r*4*u**2)/(1+r); Tp=m*v2beam/(3*k)
f=sp.factor(sp.simplify(Tp/T)); print('ratio',f)
print('f-1 =',sp.factor(f-1))
cr=sp.solve(sp.diff(f,r),r); print('crit',cr)
for c in cr:
    print('fmax',sp.nsimplify(sp.simplify(f.subs(r,c))), sp.N(f.subs(r,c)), 'x=',sp.simplify(c/2))
# Monte Carlo check that effusion flux ∝ n v for isotropic speeds: simulate particles in a box hitting a hole region
rng=np.random.default_rng(1)
def flux(v,N=400000,L=1.0,t=1.0):
    pos=rng.random((N,3))*L
    d=rng.normal(size=(N,3)); d/=np.linalg.norm(d,axis=1)[:,None]
    vx=v*d[:,0]
    # count crossings of wall x=L in time t (unfolded, ignoring hole size, wall area L^2) -> per unit area
    # use time t small: crossings = those with vx>0 and L-x < vx t
    tt=0.05/v
    return np.sum((vx>0)&(L-pos[:,0]<vx*tt))/tt
f1=flux(1.0); f2=flux(2.0); print('MC flux ratio v=2 vs v=1:',f2/f1,' theory 2; per-unit n v/4 ->',f1/400000, 0.25)
# Maxwell check: beam <v^2> = 4kT/m
v=sp.symbols('v',positive=True); a=sp.symbols('a',positive=True)
fM=v**2*sp.exp(-v**2/a**2)
gas=sp.integrate(v**2*fM,(v,0,sp.oo))/sp.integrate(fM,(v,0,sp.oo))
beam=sp.integrate(v**3*fM,(v,0,sp.oo))/sp.integrate(v*fM,(v,0,sp.oo))
print('Maxwell beam/gas',sp.simplify(beam/gas))
# selector: transit time l/v, disk rotates omega*l/v = phi + 2 pi n
l,w,ph=sp.symbols('l omega phi',positive=True); n=sp.symbols('n',integer=True,nonnegative=True)
print('v(phi)=',sp.solve(sp.Eq(w*l/v,ph+2*sp.pi*n),v))
