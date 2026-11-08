# N18-1: stepwise / ramp charging of a capacitor, Joule heat and source work
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp, quad
C,V,R,N,T,t=sp.symbols('C V R N T t',positive=True)
# step from (j-1)V/N to jV/N : heat = 1/2 C (V/N)^2 (energy balance), source work = (jV/N)*C V/N
j=sp.symbols('j',integer=True,positive=True)
H=sp.summation(C*(V/N)**2/2,(j,1,N)); W=sp.summation((j*V/N)*C*V/N,(j,1,N))
print('heat',sp.simplify(H),'work',sp.factor(W), 'stored', sp.simplify(W-H))
# discharge in N steps: heat again CV^2/(2N); energy returned to source:
Wret=sp.summation(((j-1)*V/N)*C*V/N,(j,1,N)); print('returned',sp.factor(Wret), 'check', sp.simplify(C*V**2/2 - H - Wret))
# numeric check of step charging with ODE (R=C=V=1, N=4)
def heat_profile(src,Tend):
    f=lambda tt,q:[(src(tt)-q[0])/1.0]
    s=solve_ivp(f,[0,Tend],[0.0],max_step=0.005,rtol=1e-10,atol=1e-12,dense_output=True)
    I=lambda tt:(src(tt)-s.sol(tt)[0])
    return quad(lambda tt:I(tt)**2,0,Tend,limit=2000,points=None)[0]
Nn=4; src=lambda tt: min(int(tt//30)+1,Nn)/Nn
print('step N=4 numeric heat',heat_profile(src,30*Nn+30),'theory',1/(2*Nn))
# ramp 0->V over T then hold. exact heat:
tau=R*C
I1=C*V/T*(1-sp.exp(-t/tau))
H1=sp.integrate(R*I1**2,(t,0,T))
IT=I1.subs(t,T); H2=R*IT**2*tau/2   # after T: exponential decay, heat = 1/2 R I^2 tau
Hr=sp.simplify(H1+H2)
print('ramp heat exact',Hr)
print('leading for T>>RC', sp.limit(Hr*T/(R*C**2*V**2),T,sp.oo))
for Tn in [10,50,200]:
    print(Tn, heat_profile(lambda tt:min(tt/Tn,1.0),Tn+40), float(Hr.subs({R:1,C:1,V:1,T:Tn})), 1/Tn)
