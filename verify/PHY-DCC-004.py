# N29-1: two cells in parallel feeding load R
import sympy as sp
E1,E2,r1,r2,R,I=sp.symbols('E1 E2 r1 r2 R I',positive=True)
V=sp.symbols('V')
# terminal voltage V, load current I: (E1-V)/r1+(E2-V)/r2 = I
Vsol=sp.solve((E1-V)/r1+(E2-V)/r2-I,V)[0]; print('V(I)=',sp.collect(sp.expand(Vsol),I))
E0=(E1*r2+E2*r1)/(r1+r2); r0=r1*r2/(r1+r2)
print(sp.simplify(Vsol-(E0-r0*I)))
# cell 2 charged when V>E2 ; V = E0 R/(R+r0)
Rc=sp.solve(sp.Eq(E0*R/(R+r0),E2),R)[0]; print('Rc',sp.simplify(Rc))
Pmax=E0**2/(4*r0); P1=E1**2/(4*r1)
cond=sp.factor(sp.simplify((Pmax-P1)*4*r1*r2*(r1+r2))); print('Pmax-P1 sign ~',cond)
k=sp.symbols('k',positive=True)  # k=E2/E1
print(sp.solve(sp.Eq(r1*k**2+2*r2*k-r2,0),k))
# numeric case: E1=1.6, r1=0.4, E2=1.2? r2=0.8
val={E1:sp.Rational(16,10),r1:sp.Rational(4,10),E2:sp.Rational(9,10),r2:sp.Rational(4,10)}
print('E0,r0',[sp.nsimplify(z.subs(val)) for z in (E0,r0)],'Rc',sp.nsimplify(Rc.subs(val)),'Pmax',sp.nsimplify(Pmax.subs(val)),'P1',sp.nsimplify(P1.subs(val)))
# at max power R=r0: current in cell2
Rm=r0.subs(val); Vm=(E0*Rm/(Rm+r0)).subs(val); print('V at Pmax',Vm,'I2=(E2-V)/r2',((E2-Vm)/r2).subs(val))
print('threshold k for r1=r2',sp.N(sp.sqrt(2)-1))
