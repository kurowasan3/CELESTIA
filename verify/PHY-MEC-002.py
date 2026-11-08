# T20111-1: rod-connected A (floor, rough) and B released from top; slip onset
import sympy as sp, numpy as np
th,m,M,g,l,mu=sp.symbols('theta m M g l mu',positive=True)
c=sp.cos(th); s=sp.sin(th)
v2=2*g*l*(1-c)
# radial eq for B (A at rest): toward A (center) direction: m v^2/l = m g cos(th) - F  (F: push from rod, outward)
F=sp.simplify(m*g*c - m*v2/l); print('F =',sp.factor(F))
# forces on A from rod (compression F>0 pushes A along B->A direction): horizontal -F sin (left), vertical -F cos (down)
H=-F*s; N=M*g+F*c
ratio=sp.simplify(H/N)
k=sp.symbols('k',positive=True)
r=sp.simplify(ratio.subs(M,k*m))
print('H/N =',r)
# (3) slip left at cos=4/5: mu = |H|/N
cc=sp.Rational(4,5); ss=sp.Rational(3,5)
mu3=sp.simplify((-r).subs({sp.cos(th):cc,sp.sin(th):ss})); print('mu(c=4/5)=',sp.factor(mu3))
x=sp.symbols('x',positive=True)
hfun=(3*x-2)*sp.sqrt(1-x**2)/(k+(3*x-2)*x)
cond=sp.factor(sp.together(sp.diff(hfun,x).subs(x,cc))); print('dh/dc at 4/5:',cond)
# numeric scan for k=1/5
kk=sp.Rational(1,5); hf=sp.lambdify(x,hfun.subs(k,kk))
xs=np.linspace(0.8,1,4001)[1:]; print('max h on (4/5,1):',hf(xs).max(),' mu=',float(mu3.subs(k,kk)), '=6/13',6/13)
# (4) liftoff: N=0 in tension phase: k=(2-3c)c
print('liftoff c:',sp.solve(sp.Eq(k,(2-3*x)*x),x))
# for k=1/5: compute pull/N as c decreases from 2/3 -> diverges before liftoff
gf=sp.lambdify(x,((2-3*x)*sp.sqrt(1-x**2)/(kk-(2-3*x)*x)))
cl=float((1+sp.sqrt(1-3*kk))/3); xs2=np.linspace(2/3-1e-6,cl+1e-6,5); print('c_lift',cl,'ratio',gf(xs2))
