# N22-1: charge pump with C1, C2 (general), and a "load" capacitor C3 that takes charge each cycle
import sympy as sp
from fractions import Fraction as Fr
C,V0,c1,c2,x=sp.symbols('C V0 c1 c2 x',positive=True)
# recurrence (no load)
xn=sp.symbols('x_n')
xnext=sp.solve(sp.Eq(c1*(x-V0)+c2*x, c1*V0+c2*xn),x)[0]; print('x_{n+1}=',sp.simplify(xnext))
print('fixed point', sp.solve(sp.Eq(xnext.subs(xn,x),x),x))
# simulate exact cycle with fractions
def sim(C1,C2,C3,ncyc,V=Fr(1)):
    Q1=Fr(0);Q2=Fr(0);Wb=Fr(0);Hpump=Fr(0);Hload=Fr(0);Eload=Fr(0);hist=[]
    E=lambda Cc,Q: Q*Q/(2*Cc)
    for n in range(ncyc):
        # left: C1 connected across battery
        dq=C1*V-Q1; Wb+=V*dq; Hpump+= V*dq-(E(C1,C1*V)-E(C1,Q1)); Q1=C1*V
        # right: loop battery + C1 -> C2 ; V0+Q1'/C1 = Q2'/C2 ; Q1'+Q2'=Q1+Q2
        S=Q1+Q2
        Q2n=(V+S/C1)/(1/C1+1/C2); Q1n=S-Q2n
        dq=Q2n-Q2; Wb+=V*dq; Hpump+=V*dq-(E(C1,Q1n)+E(C2,Q2n)-E(C1,Q1)-E(C2,Q2))
        Q1,Q2=Q1n,Q2n; hist.append(Q2/C2)
        if C3:
            Etot=E(C2,Q2); Q2n=Q2*C2/(C2+C3); Q3=Q2*C3/(C2+C3)
            Hload+=Etot-E(C2,Q2n)-E(C3,Q3); Eload+=E(C3,Q3); Q2=Q2n
    return hist,Wb,Hpump,Hload,Eload
h,*_=sim(Fr(1),Fr(1),None,4); print('C1=C2 V_AB:',[str(v) for v in h])
h,*_=sim(Fr(1),Fr(3),None,14); print('C2=3C1:',[float(v) for v in h])
for n,v in enumerate(h,1):
    if v>Fr(19,10): print('first n exceeding 1.9V0:',n); break
print('formula 2(1-(3/4)^n):',[2*(1-0.75**n) for n in (10,11)])
# load case C1=C2=C3=1
h,Wb,Hp,Hl,El=sim(Fr(1),Fr(1),Fr(1),200)
print('steady V_AB just after pump',float(h[-1]))
# per-cycle in steady state: difference of 2 long runs
h2,Wb2,Hp2,Hl2,El2=sim(Fr(1),Fr(1),Fr(1),201)
print('per cycle: Wb',float(Wb2-Wb),'Hpump',float(Hp2-Hp),'Hload(share)',float(Hl2-Hl),'Eload',float(El2-El))
print('expected 4/3, 2/3, 4/9, 2/9')
# symbolic steady state with C3
xs=sp.solve(sp.Eq(x,(2*c1*V0+c2*x*c2/(c2+C))/(c1+c2)),x)[0]; print('steady general', sp.factor(xs))
