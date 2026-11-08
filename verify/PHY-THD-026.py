# K20133-1 verification
import sympy as sp
p,p1,V1=sp.symbols('p p1 V1',positive=True)
a=sp.Rational(5,3); b=sp.Rational(7,5)
VA=V1*(p1/p)**(1/a); VB=V1*(p1/p)**(1/b)
print('(1) VA=',VA,' VB=',VB, ' VB/VA=',sp.simplify(VB/VA))
ps=sp.solve(sp.Eq(VB,2*V1),p)[0]; print('(2) ps=',ps,float(ps/p1),' VA(ps)=',sp.simplify(VA.subs(p,ps)),float(VA.subs(p,ps)/V1))
pf=p1/32
VAf=sp.simplify(VA.subs(p,pf)); print('(3) VAf=',VAf)
dUA=sp.Rational(3,2)*(pf*VAf-p1*V1); dUB=sp.Rational(5,2)*(ps*2*V1-p1*V1)
Wtot=sp.simplify(-(dUA+dUB)); print('    W on piston =',Wtot, float(Wtot/(p1*V1)))
# independent check: integrate p dV_total numerically with stopper
import numpy as np
from scipy.integrate import quad
p1n=1.0;V1n=1.0
psn=float(ps.subs(p1,1))
# work = ∫ p dV_total ; parametrize by p. phase1: dVtot = dVA+dVB; phase2: only VA
dVA=lambda pp: -(3/5)*V1n*p1n**(3/5)*pp**(-8/5)
dVB=lambda pp: -(5/7)*V1n*p1n**(5/7)*pp**(-12/7)
w1=quad(lambda pp: -pp*(dVA(pp)+dVB(pp)), psn, 1)[0]
w2=quad(lambda pp: -pp*dVA(pp), 1/32, psn)[0]
print('    numeric ∫p dV =',w1+w2)
# (4) equal work without stopper
r=sp.Rational(1,32)
WA=lambda VA1: VA1*p1/(a-1)*(1-r**((a-1)/a))
WB=lambda VB1,g: VB1*p1/(g-1)*(1-r**((g-1)/g))
x=sp.symbols('x',positive=True)  # VA1/VB1
sol=sp.solve(sp.Eq(WA(x),WB(1,b)),x)[0]; print('(4) VA1/VB1=',sp.nsimplify(sol),float(sol),' check form',float(sp.Rational(20,9)*(1-2**sp.Rational(-10,7))))
print('    WA coef',sp.simplify(WA(1)/p1),' WB coef',float(WB(1,b)/p1))
# (5) isothermal B
WBiso=p1*sp.log(32)
sol2=sp.solve(sp.Eq(WA(x),WBiso),x)[0]; print('(5) VA1/VB1 iso =',sp.simplify(sol2),float(sol2))
g=sp.symbols('g',positive=True)
print('    limit g->1 of WB:',sp.limit(WB(1,g),g,1))
# monotonic in gamma check
print('    WB(gamma) for 1.01,1.4,1.67:',[float(WB(1,sp.nsimplify(gg))/p1) for gg in [1.01,1.4,5/3]])
