# Q059-1: hovering power from momentum flux, scaling
import sympy as sp
M,g,rho,A,w,k,beta=sp.symbols('M g rho A w k beta',positive=True)
wsol=sp.solve(sp.Eq(rho*A*w*w,M*g),w)[0]
P=sp.simplify(sp.Rational(1,2)*rho*A*w*wsol**2).subs(w,wsol)
P=sp.simplify(P); print('w',wsol,' P',P)
# alternative: P = (1/2) thrust * w
print('check P=Mg w/2:',sp.simplify(P-M*g*wsol/2))
Pk=P.subs({M:M*k**3,A:A*k**2}); r=sp.simplify(Pk/P); print('ratio',r, ' k=1/2:',sp.N(r.subs(k,sp.Rational(1,2))), 1/sp.N(r.subs(k,sp.Rational(1,2))))
kmax=sp.solve(sp.Eq(beta*k**3,r),k); print('kmax',kmax, [sp.N(v.subs(beta,1.5)) for v in kmax])
# altitude: rho'=0.70 rho, max mass M' with power beta*P
Mp=sp.symbols('Mp',positive=True)
eq=sp.Eq(P.subs({M:Mp,rho:sp.Rational(7,10)*rho}),beta*P)
s=sp.solve(eq,Mp)[0]; print('Mp/M',sp.simplify(s/M), sp.N((s/M).subs(beta,1.5)))
