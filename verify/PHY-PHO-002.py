# N56-1 photoelectric charging of an isolated parallel-plate capacitor
import sympy as sp
h,nu,W,e,C,eps0,S,d,x=sp.symbols('h nu W e C epsilon_0 S d x',positive=True)
V1=(h*nu-W)/e
Q1=C*V1; N1=sp.simplify(Q1/e); print('N1=',N1)
# halve gap slowly with light on: V stays V1 (emission tops up), C->2C
Q2=2*C*V1; print('extra electrons',sp.simplify((Q2-Q1)/e))
# external work: F = Q^2/(2 eps0 S), Q=eps0 S V1/x ; agent force outward, displacement inward
Fx=(eps0*S*V1/x)**2/(2*eps0*S)
Wext=-sp.integrate(Fx,(x,d/2,d)); Wext=sp.simplify(Wext.subs(S,C*d/eps0)); print('W_ext=',Wext)
# energy balance: dU = Wext + (charge moved)*V1
dU=sp.Rational(1,2)*2*C*V1**2-sp.Rational(1,2)*C*V1**2
print('balance',sp.simplify(dU-(Wext+(Q2-Q1)*V1)))
# back to d in dark: Q fixed 2CV1 -> V=2V1 ; threshold nu2 for renewed flow
nu2=sp.solve(sp.Eq(h*x-W,2*e*V1),x)[0]; print('nu2=',sp.simplify(nu2))
# with nu2 light, after long time voltage -> (h nu2' - W)/e; check with nu3 = 3nu-?; also Wext for return d/2 -> d at fixed Q
Q=Q2; Wret=sp.integrate(Q**2/(2*eps0*S),(x,d/2,d)); print('W return',sp.simplify(Wret.subs(S,C*d/eps0)))
