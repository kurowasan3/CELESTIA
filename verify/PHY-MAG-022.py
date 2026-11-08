# Q172-1 thick solenoid winding: heat independent of wire size; two-layer optimum; continuous limit
import sympy as sp
a,b,c,d,d1,d2,l,rho,mu,B,I=sp.symbols('a b c d d1 d2 ell rho mu0 B I',positive=True)
# single wire, square cross-section d x d
N=(b-a)*l/d**2; Lw=N*sp.pi*(a+b); Rw=rho*Lw/d**2
Ireq=B*l/(mu*N); P=sp.simplify(Ireq**2*Rw); V=sp.simplify(Ireq*Rw)
print('P=',sp.factor(P),' V=',sp.factor(V),' I=',sp.factor(Ireq))
# two layers series, same I
N1=(c-a)*l/d1**2; N2=(b-c)*l/d2**2
R1=rho*N1*sp.pi*(a+c)/d1**2; R2=rho*N2*sp.pi*(c+b)/d2**2
Ib=B*l/(mu*(N1+N2)); P2=sp.simplify(Ib**2*(R1+R2))
x=sp.symbols('x',positive=True)  # x = d1^2/d2^2
P2x=sp.simplify(P2.subs(d1,sp.sqrt(x)*d2))
xs=sp.solve(sp.diff(P2x,x),x); print('x opt=',[sp.simplify(s) for s in xs])
Pmin=sp.simplify(P2x.subs(x,xs[0])); print('Pmin=',sp.factor(Pmin))
K=sp.pi*rho*l*(B/mu)**2
w1=(c+a)/(c-a); w2=(b+c)/(b-c)
print('check formula', sp.simplify(Pmin-K/(1/w1+1/w2)))
sub={c:2*a,b:3*a}
print('uniform/K',sp.simplify((P/K).subs(sub)),' opt/K',sp.simplify((Pmin/K).subs(sub)),' ratio',sp.simplify((Pmin/P).subs(sub)), ' d1/d2=',sp.sqrt(sp.simplify(xs[0].subs(sub))))
# continuous limit: J(r)=A/r ; B/mu=int J dr ; P/(rho l)= int J^2 2 pi r dr
r,A=sp.symbols('r A',positive=True)
Aval=sp.solve(sp.Eq(sp.integrate(A/r,(r,a,b)),B/mu),A)[0]
Pc=sp.simplify(rho*l*sp.integrate((Aval/r)**2*2*sp.pi*r,(r,a,b)))
print('Pcont/K=',sp.simplify(Pc/K), float((Pc/K).subs(b,3*a)))
# check with N-layer discrete optimum (generic) numerically approaches
import numpy as np
for n in [2,10,1000]:
    edges=np.linspace(1,3,n+1); w=(edges[1:]+edges[:-1])/(edges[1:]-edges[:-1]); print(n,1/np.sum(1/w))
