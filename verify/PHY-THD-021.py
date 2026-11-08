# T20203-1: Brayton-type refrigerator cooling a finite body Y (1 mol monatomic, const volume)
import sympy as sp
TH,b,R,TY=sp.symbols('T_H b R T_Y',positive=True)
Cp=sp.Rational(5,2)*R; Cv=sp.Rational(3,2)*R
# adiabat pH -> pL=pH/b^5 : T ∝ p^{2/5}
T1=TH*(1/b**5)**sp.Rational(2,5); print('T after expansion:',sp.simplify(T1), ' V ratio:',sp.simplify((TH/b**2)/(TH)*b**5))
Tp=sp.symbols('Tp'); Tprime=sp.solve(sp.Eq(Cp*(Tp-TH/b**2),Cv*(TY-Tp)),Tp)[0]; print("T'=",sp.simplify(Tprime))
QY=sp.simplify(Cp*(Tprime-TH/b**2)); QZ=sp.simplify(Cp*(b**2*Tprime-TH)); print('QZ/QY=',sp.simplify(QZ/QY))
# work done on gas per cycle (gas returns to initial state)
W=sp.simplify(QZ-QY); print('W/QY=',sp.simplify(W/QY))
# recurrence and limit
n=sp.symbols('n',integer=True,nonnegative=True)
Tinf=TH/b**2
Tn=Tinf+(TH-Tinf)*sp.Rational(3,8)**n
print('check recurrence:',sp.simplify(Tprime.subs(TY,Tn)-Tn.subs(n,n+1)))
print('T1 =',sp.simplify(Tn.subs(n,1)))
Wtot=sp.simplify((b**2-1)*Cv*(TH-Tinf)); print('Wtot=',sp.factor(Wtot))
# brute-force numeric: iterate cycles and sum works directly
import numpy as np
THv,bv,Rv=300.0,np.sqrt(2),8.314
TYv=THv; Wsum=0; QYs=0
for i in range(200):
    T1v=THv/bv**2
    Tpv=(5*T1v+3*TYv)/8
    qy=2.5*Rv*(Tpv-T1v)
    # adiabatic compression work on gas = Cv*(b^2 Tp - Tp); expansion work by gas = Cv*(TH - T1); isobaric works
    w_on=1.5*Rv*(bv**2*Tpv-Tpv) - 1.5*Rv*(THv-T1v) - Rv*(Tpv-T1v) + Rv*(bv**2*Tpv-THv)
    Wsum+=w_on; TYv=Tpv
print('numeric limit T_Y',TYv,' pred',THv/bv**2,' Wtot numeric',Wsum,' pred',float(Wtot.subs({TH:THv,b:bv,R:Rv})))
print('b^2=2 :',sp.simplify(Wtot.subs(b,sp.sqrt(2))))
# specific case p_H/p_L = 32  ->  b^5=32, b^2=4
sub={b:2}
print('T after expansion:',sp.simplify(T1.subs(sub)),' T_Y1:',sp.simplify(Tn.subs(n,1).subs(sub)),' Tinf:',Tinf.subs(sub))
print('W_total:',sp.simplify(Wtot.subs(sub)),' total heat from Y:',sp.simplify((Cv*(TH-Tinf)).subs(sub)),' heat to Z total:',sp.simplify((b**2*Cv*(TH-Tinf)).subs(sub)))
print('first-cycle W:',sp.simplify((W.subs(TY,TH)).subs(sub)))
import math
for N in range(1,10):
    if (3/8)**N*(3/4) <= 0.01*(1/4): print('cycles needed for within 1% of T_inf:',N); break
# V ratios (pV=RT): state0 V0=R TH/pH ; state1 = 8 V0
print('V1/V0=',sp.simplify((TH/b**2)/(TH)*b**5).subs(sub))
