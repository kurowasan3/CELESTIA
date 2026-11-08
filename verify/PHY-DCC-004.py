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
# old numeric example (kept as a consistency check of (2)-(4))
val={E1:sp.Rational(16,10),r1:sp.Rational(4,10),E2:sp.Rational(9,10),r2:sp.Rational(4,10)}
print('E0,r0',[sp.nsimplify(z.subs(val)) for z in (E0,r0)],'Rc',sp.nsimplify(Rc.subs(val)),'Pmax',sp.nsimplify(Pmax.subs(val)),'P1',sp.nsimplify(P1.subs(val)))
# at max power R=r0: current in cell2
Rm=r0.subs(val); Vm=(E0*Rm/(Rm+r0)).subs(val); print('V at Pmax',Vm,'I2=(E2-V)/r2',((E2-Vm)/r2).subs(val))
print('threshold k for r1=r2',sp.N(sp.sqrt(2)-1))

# ---- (5) rewritten ----
Vmax=E0/2
cond=sp.simplify((Vmax-E2)*2*(r1+r2))
print('(5a) V-E2 at Pmax, times 2(r1+r2):',sp.expand(cond),' -> charged iff E1 r2 > E2 (r1+2 r2)')
assert sp.expand(cond-(E1*r2-E2*(r1+2*r2)))==0
# consistent with (2): r0 > r1 E2/(E1-E2)
assert sp.simplify(sp.factor((r0-Rc)*(E1-E2)*(r1+r2))-(r1*(E1*r2-E2*(r1+2*r2))))==0
# (5b) random check: whenever Pmax>PA, cell B discharges at R=r0 (I_B>0); also A's external power <= PA always
import random
cnt=0
for _ in range(200000):
    e1=random.uniform(0.1,5); e2=random.uniform(0.01,e1*0.999); a=random.uniform(0.01,5); b=random.uniform(0.01,5)
    E0n=(e1*b+e2*a)/(a+b); r0n=a*b/(a+b); V=E0n/2
    IB=(e2-V)/b; IA=(e1-V)/a
    assert V*IA<=e1**2/(4*a)+1e-12
    if E0n**2/(4*r0n)>e1**2/(4*a):
        cnt+=1; assert IB>0
print('(5b) cases with Pmax>PA:',cnt,' all have B discharging')
k_=sp.symbols('k',positive=True)
print('     r1=r2: charge iff k<1/3 ; (4) needs k>sqrt2-1 =',sp.N(sp.sqrt(2)-1),' > 1/3 :',bool(sp.sqrt(2)-1>sp.Rational(1,3)))
print('DCC-004 (5) OK')
