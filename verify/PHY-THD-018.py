# T20043-1 connected piston (areas 2S0 | S0), gas2 under weighted piston B (area S), adiabatic, then stopper
import sympy as sp
from scipy.integrate import solve_ivp
import numpy as np
R,T0,M,g,S,S0,h,x,u=sp.symbols('R T0 M g S S0 h x u',positive=True)
P2=M*g/S; P1=P2*S0/(2*S0)   # P1*2S0 = P2*S0
V1=R*T0/P1; V2=R*T0/P2
print('V1,V2=',sp.simplify(V1),sp.simplify(V2))
# stage1: gas2 adiabatic quasi-static at const P: (3/2)P dV2 = -P dV2 -> dV2=0
dV2=sp.symbols('dV2'); print('dV2 sol:',sp.solve(sp.Eq(sp.Rational(3,2)*P2*dV2,-P2*dV2),dV2))
xs=sp.solve(sp.Eq(S0*x,S*h),x)[0]; print('x=',xs)
dV1=2*S0*xs; T1=sp.simplify(P1*(V1+dV1)/R); print('T1=',T1)
Q1=sp.simplify(sp.Rational(3,2)*R*(T1-T0)+P1*dV1); print('Q1=',Q1)
# stage2: gas2 adiabatic from V2 to V2/8 (gamma=5/3) -> T2=4T0, P2'=32P2
g53=sp.Rational(5,3)
P2p=P2*8**g53; T2p=T0*8**(g53-1); print('P2p/P2,T2p=',sp.simplify(P2p/P2),T2p)
P1p=P2p/2
V1p=V1+dV1+2*(V2-V2/8)
T1p=sp.simplify(P1p*V1p/R); print('T1p=',sp.expand(T1p))
Q2=sp.simplify(sp.Rational(3,2)*R*(T1p-T1)+sp.Rational(3,2)*R*(T2p-T0)); print('Q2=',sp.factor(sp.expand(Q2)))
# Q2 via work integral: W1 = int P1 dV1 = int P2 dV2(-) = Delta U2
W1=sp.Rational(3,2)*R*(T2p-T0)
Q2b=sp.Rational(3,2)*R*(T1p-T1)+W1
print('check same:',sp.simplify(Q2-Q2b)==0)
# numeric: integrate stage 2 quasi-statically with u = S0*x (volume swept into container2)
vals={R:8.31,T0:300.0,M:10.0,g:9.8,S:0.01,S0:0.02,h:0.05}
P2n=float(P2.subs(vals)); V2n=float(V2.subs(vals)); V1s=float((V1+dV1).subs(vals))
def P2u(uu): return P2n*(V2n/(V2n-uu))**(5/3)
us=np.linspace(0,7*V2n/8,200001)
P1s=P2u(us)/2
W=np.trapezoid(P1s*2,us)   # P1 dV1 with dV1 = 2 du
T1n=P1s[-1]*(V1s+2*us[-1])/8.31
Q2n=1.5*8.31*(T1n-float(T1.subs(vals)))+W
print('Q2 numeric',Q2n,' formula',float(Q2.subs(vals)))
assert abs(Q2n/float(Q2.subs(vals))-1)<1e-6
# check gas2 energy: W equals Delta U2
print('W1 numeric',W,' DeltaU2',1.5*8.31*3*300)
print('OK')
