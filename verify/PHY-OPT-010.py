# N09-1: fringe count measures density, not temperature
import sympy as sp
T0,T,w,lam,d0,g=sp.symbols('T0 T w lambda delta0 gamma',positive=True)
# n-1 = K*rho ; constant pressure: rho ∝ 1/T
N=(d0 - d0*T0/T)*w/lam; print("N(T)=",sp.simplify(N), " limit:",sp.limit(N,T,sp.oo))
vals={d0:2.53e-4,w:0.100,lam:633e-9,T0:300}
C=(d0*w/lam).subs(vals); print("delta0*w/lam =",C)
Tsol=sp.solve(sp.Eq(N.subs(vals),10),T); print("T =",Tsol)
R=8.31; nu=41.6/(2.5*R*(Tsol[0]-300)); print("nu =",nu)
# alternative: first law: Q = dU + p dV = 3/2 nu R dT + nu R dT
print("Q check",(1.5+1)*0.0200*R*(Tsol[0]-300))
# adiabatic halving: rho doubles -> N = d0*w/lam*(2-1)
print("adiabatic N =",C, " T =",300*2**(sp.Rational(2,3)), float(300*2**(sp.Rational(2,3))))
