# O20173-2 sand grains lifted onto piston: hand work = ∫V dp ; isothermal vs adiabatic
import sympy as sp, numpy as np
from scipy.integrate import quad
p0,V0,p1,p=sp.symbols('p0 V0 p1 p',positive=True)
x=sp.symbols('x',positive=True) # p1/p0
# isothermal
Viso=p0*V0/p
Wiso=sp.integrate(Viso,(p,p0,p1)); print("W_hand iso=",sp.simplify(Wiso))
V1i=p0*V0/p1; Esand_i=(p1-p0)*V1i
print("E_sand iso=",sp.simplify(Esand_i))
Qout_i=sp.simplify(Wiso+p0*(V0-V1i)-Esand_i-0); print("Q_out iso =",Qout_i)
# adiabatic gamma=5/3
Vad=V0*(p0/p)**sp.Rational(3,5)
Wad=sp.simplify(sp.integrate(Vad,(p,p0,p1))); V1a=Vad.subs(p,p1)
print("W_hand ad=",Wad," vs 5/2(p1V1-p0V0)=",sp.simplify(Wad-sp.Rational(5,2)*(p1*V1a-p0*V0)))
dU=sp.Rational(3,2)*(p1*V1a-p0*V0); Esand_a=(p1-p0)*V1a; Watm=p0*(V0-V1a)
print("energy balance ad: W_hand + W_atm - E_sand - dU =",sp.simplify(Wad+Watm-Esand_a-dU))
Won=-sp.integrate(p0*V0**sp.Rational(5,3)*sp.Symbol('V')**sp.Rational(-5,3),(sp.Symbol('V'),V0,V1a))
print("W_hand/W_on_gas =",sp.simplify(Wad/Won))
# then cool at p1 to V1i; total heat released
Qcool=sp.Rational(5,2)*p1*(V1a-V1i)
# overall energy for route2: W_hand_ad + W_atm(total, V0->V1i) - E_sand(final) = Qout_total (dU=0)
Qout_2=sp.simplify(Wad+p0*(V0-V1i)-Esand_i)
print("Qout route2 via balance =",Qout_2, " check vs (work on gas during ad? ) -> Qcool=",sp.simplify(Qcool))
print("difference Qout2-Qout_iso =",sp.simplify(Qout_2-Qout_i - (Wad-Wiso)))
# check Qout2 == Qcool  (adiabatic part no heat)
print("Qout2 - Qcool =",sp.simplify(Qout_2-Qcool))
# hmm need: sand stays on piston during cooling, piston descends -> sand loses PE, atmosphere works
# numbers x=2
for name,expr in [("Wiso",Wiso),("Wad",Wad),("Esand_i",Esand_i),("Qcool",Qcool)]:
    print(name, sp.N(sp.expand_log(sp.simplify(expr.subs(p1,2*p0)/(p0*V0)),force=True)))
print("2^(2/5)=",2**0.4)
# numeric: discrete grains simulation, adiabatic
N=200000; pp=1.0; VV=1.0; Wh=0.0; dp=1.0/N
for i in range(N):
    Wh+=VV*dp   # lift grain of mass S dp/g to height V/S
    pp+=dp; VV=(1.0/pp)**0.6
print("discrete adiabatic hand work",Wh)
