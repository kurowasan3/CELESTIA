# K20253-1: gas with U = alpha P V and u=u(T) -> P ∝ T^(alpha+1)
import sympy as sp, numpy as np
from scipy.integrate import quad
a,P,V,T,c,dT,V1,V2=sp.symbols('alpha P V T c dT V1 V2',positive=True)
# (2) adiabat: dU = -P dV with U = a P V
Pf=sp.Function('P')
sol=sp.dsolve(sp.Eq(a*(Pf(V)+V*Pf(V).diff(V)),-Pf(V)),Pf(V)); print('(2) adiabat:',sol,' -> gamma=(a+1)/a')
# (1) isothermal: P const, Q = dU + P dV = a P dV + P dV
print('(1) Q=(a+1)P(V2-V1)')
# (3) thin cycle: efficiency dP(V2-V1)/((a+1)P(V2-V1)) = dT/T  -> dP/dT=(a+1)P/T
Pt=sp.Function('P'); s=sp.dsolve(sp.Eq(Pt(T).diff(T),(a+1)*Pt(T)/T),Pt(T)); print('(3)',s)
# numerical check: finite Carnot cycle with P=c T^(a+1) has efficiency exactly 1-TL/TH; other exponent fails
def carnot_eff(alpha,x,TH=2.0,TL=1.0,VA=1.0,VB=3.0):
    Pof=lambda T: T**x
    g=(alpha+1)/alpha
    PH,PL=Pof(TH),Pof(TL)
    # isothermal at TH from VA to VB ; adiabat from VB to P=PL ; isothermal back ; adiabat to VA
    VC=VB*(PH/PL)**(1/g); VD=VA*(PH/PL)**(1/g)
    QH=(alpha+1)*PH*(VB-VA); QL=(alpha+1)*PL*(VC-VD)
    Wad1=quad(lambda v: PH*(VB/v)**g,VB,VC)[0]; Wad2=quad(lambda v: PL*(VD/v)**g,VD,VA)[0]
    W=PH*(VB-VA)+Wad1-PL*(VC-VD)+Wad2
    return W/QH, QH/TH-QL/TL
for al in (1,2,3,1.5):
    print(f' alpha={al}: x=alpha+1 -> eff={carnot_eff(al,al+1)[0]:.6f} (Carnot 0.5), Clausius sum={carnot_eff(al,al+1)[1]:.2e};  x=alpha -> eff={carnot_eff(al,al)[0]:.6f}')
# (4) free expansion: U const -> u halves -> T'=T*2^(-1/(a+1)); slow adiabatic: T V^(1/a) const
x=a+1
Tad=sp.simplify((2)**(-1/a)); Tfree=2**(-1/(a+1))
print('(4) free:',Tfree,' adiabatic:',Tad)
# check adiabatic T V^(1/a): P V^((a+1)/a) const and P ∝ T^(a+1)
print('   T^(a+1) V^((a+1)/a) = (T V^(1/a))^(a+1)  ok')
for al in (1,2,3): print('   alpha',al,'free',2**(-1/(al+1)),'adiab',2**(-1/al))
# (5) spring: P = beta V, heat from T0 to 2T0: Q = a Δ(PV) + ∫P dV = (a+1/2)(P2V2-P1V1); P2V2=4^(a+1) P1V1
b,Vv=sp.symbols('beta Vv',positive=True)
V1s=sp.symbols('V_1',positive=True); k=2**(a+1)
Q=a*(b*(k*V1s)**2-b*V1s**2)+sp.integrate(b*Vv,(Vv,V1s,k*V1s))
print('(5) Q/(P1V1)=',sp.simplify(Q/(b*V1s**2)),' expected (a+1/2)(4^(a+1)-1):',sp.simplify(Q/(b*V1s**2)-(a+sp.Rational(1,2))*(4**(a+1)-1))==0)
# check T doubling consistent: P=βV and P=cT^(a+1): V∝T^(a+1) -> V2/V1 = 2^(a+1)
print('   alpha=3:',float(((3+0.5)*(4**4-1))))
