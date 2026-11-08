# Q062-1: two connected soap bubbles, isothermal, determine gamma; sensitivity
import sympy as sp
p0,g,R1,R2,R3=sp.symbols('p0 gamma R1 R2 R3',positive=True)
pV=lambda R:(p0+4*g/R)*R**3
gam=sp.solve(sp.Eq(pV(R1)+pV(R2),pV(R3)),g)[0]; print('gamma =',sp.factor(gam))
num={p0:1.0e5,R1:1.0e-2,R2:2.0e-2}
gval=2.5e-2
import mpmath as mp; mp.mp.dps=40
R3s=sp.nsolve((pV(R1)+pV(R2)-pV(R3)).subs({p0:sp.Integer(10)**5,R1:sp.Rational(1,100),R2:sp.Rational(2,100),g:sp.Rational(25,1000)}),R3,sp.Float('0.0208',40),prec=40)
Rstar=sp.N((sp.Rational(9,10**6))**sp.Rational(1,3),40)
print('R3',R3s,' R*',Rstar,' delta',R3s-Rstar)
# approximate formula delta ~ eps (R1^2+R2^2-R3^2)/(3 R3^2), eps=4g/p0
eps=4*0.025/1e5
print('approx delta',eps*(1e-4+4e-4-float(Rstar)**2)/(3*float(Rstar)**2))
# sensitivity: d gamma/d R3 at the solution
dg=sp.diff(gam,R3).subs({p0:sp.Integer(10)**5,R1:sp.Rational(1,100),R2:sp.Rational(2,100),R3:R3s})
print('dgamma/dR3',sp.N(dg),'  dR3 for 10% of gamma:',sp.N(0.1*0.025/abs(dg)))
# check air flows small->large: pressure larger in small
print('excess p small,large',4*0.025/0.01,4*0.025/0.02)
