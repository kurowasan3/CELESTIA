# N27-1: two rods in series, Drude-like model (accelerate from rest for tau then stop)
import sympy as sp
e,m,S,l1,l2,n1,n2,t1,t2,V,eps=sp.symbols('e m S l1 l2 n1 n2 tau1 tau2 V epsilon0',positive=True)
E1,E2=sp.symbols('E1 E2',positive=True)
vbar=lambda E,t: e*E*t/(2*m)
sol=sp.solve([e*n1*S*vbar(E1,t1)-e*n2*S*vbar(E2,t2), E1*l1+E2*l2-V],[E1,E2],dict=True)[0]
E1s,E2s=sp.simplify(sol[E1]),sp.simplify(sol[E2]); print('E1',E1s,'E2',E2s)
I=sp.simplify(e*n1*S*vbar(E1s,t1)); print('I',I)
rho=lambda n,t: 2*m/(e**2*n*t)
Rtot=rho(n1,t1)*l1/S+rho(n2,t2)*l2/S
print('I check',sp.simplify(I-V/Rtot))
Q=sp.simplify(eps*S*(E2s-E1s)); print('junction charge',sp.factor(Q))
print('drift ratio v2/v1',sp.simplify(vbar(E2s,t2)/vbar(E1s,t1)))
# heat per unit volume by collisions: n*(1/t)*(1/2)m(eEt/m)^2 ; compare E*j
for (n,t,Ei) in [(n1,t1,E1s),(n2,t2,E2s)]:
    p_coll=n*(1/t)*sp.Rational(1,2)*m*(e*Ei*t/m)**2
    print('heat/vol consistent with E*j:',sp.simplify(p_coll-Ei*I/S))
print('ratio heat/vol 2:1',sp.simplify((E2s*I/S)/(E1s*I/S)))
# numeric example: copper-like rod1, rod2 with n2=n1/2, t2=t1/4  -> rho2=8 rho1
sub={n2:n1/2,t2:t1/4,l2:l1}
print('E2/E1',sp.simplify((E2s/E1s).subs(sub)),'v2/v1',sp.simplify((vbar(E2s,t2)/vbar(E1s,t1)).subs(sub)))
print('Q',sp.simplify(Q.subs(sub)))
