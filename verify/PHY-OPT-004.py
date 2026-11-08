# N03-1: film over one slit, apparatus in air and in liquid
import sympy as sp, numpy as np
nf,dl,d,l,nl,X1,X2,x=sp.symbols('n_f delta d l n_l X1 X2 x',positive=True)
# optical path difference (A path - B path) at screen coordinate x (upward toward B), whole space index nl, film on B
opd = nl*d*x/l - (nf-nl)*dl
x0=sp.solve(sp.Eq(opd,0),x)[0]; print("x0 =",sp.simplify(x0))
sol=sp.solve([sp.Eq(X1,(nf-1)*dl*l/d), sp.Eq(-X2,(nf-nl)*dl*l/(nl*d))],[nf,dl],dict=True)[0]
print({k:sp.simplify(v) for k,v in sol.items()})
vals={X1:sp.Rational(10,1000),X2:sp.Rational(125,100000),nl:sp.Rational(16,10),d:sp.Rational(5,10000),l:1}
print("n_f =",sol[nf].subs(vals)," delta =",sol[dl].subs(vals))
lam=600e-9
print("spacing air(mm)",lam*1/0.5e-3*1e3,"liquid(mm)",lam/(1.6*0.5e-3)*1e3)
print("shift in fringes in liquid:", (1.5-1.6)*1e-5/lam)
# brute-force: exact geometry intensity maximum location (zero order) in liquid
dd,L=0.5e-3,1.0
xs=np.linspace(-5e-3,5e-3,2000001)
A=np.sqrt(L**2+(xs+dd/2)**2); B=np.sqrt(L**2+(xs-dd/2)**2)
opd_n=1.6*(A-B)-(1.5-1.6)*1e-5
print("zero-OPD x (mm):", xs[np.argmin(abs(opd_n))]*1e3)
opd_a=(A-B)-(0.5)*1e-5
xs2=np.linspace(0,20e-3,2000001);A=np.sqrt(L**2+(xs2+dd/2)**2); B=np.sqrt(L**2+(xs2-dd/2)**2)
print("air zero-OPD x (mm):", xs2[np.argmin(abs((A-B)-0.5e-5))]*1e3)
