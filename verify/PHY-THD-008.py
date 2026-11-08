import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
p,V,T,R=sp.symbols('p V T R',positive=True)
nA=3*p*V/(R*T); nB=p*2*V/(R*2*T); n=nA+nB
# (1) U conservation
p1=sp.Rational(2,3)*sp.Rational(3,2)*(3*p*V+p*2*V)/(3*V)
T1=sp.simplify(p1*3*V/(n*R))
print('p1',sp.simplify(p1),'T1',T1)
# (2) U conservation with equal pressure -> p2
p2=sp.symbols('p2'); print('p2',sp.solve(sp.Eq(3*p*V+2*p*V, p2*V+p2*2*V),p2))
P2=sp.Rational(5,3)*p
TA=T*(P2/(3*p))**sp.Rational(2,5); nAp=P2*V/(R*TA); nBp=n-nAp; TB=P2*2*V/(R*nBp)
print('TA/T',sp.N(TA/T),'TB/T',sp.N(sp.simplify(TB/T)),'dn',sp.N((nA-nAp)*R*T/(p*V)))
# numeric: ODE in transferred moles, units p=V=R=T=1
def rhs(s,y):
    # y=[nA, TA, UB]; transfer dn=ds
    nA_,TA_,UB=y
    pA=nA_*TA_/1.0
    # adiabatic quasi-static for remaining gas: TA ~ pA^(2/5); d(nA)=-1 -> dpA/pA = ... use T V^(2/3) with specific volume v=1/nA
    # T v^(2/3)=const, v=V/nA -> dT/T = (2/3) dnA/nA
    dTA=(2/3)*TA_*(-1)/nA_
    dUB=2.5*TA_   # enthalpy inflow
    return [-1,dTA,dUB]
def ev(s,y):
    nA_,TA_,UB=y; return nA_*TA_ - UB/1.5/2.0
ev.terminal=True
sol=solve_ivp(rhs,[0,3],[3.0,1.0,1.5*1.0*2.0],events=ev,rtol=1e-10,atol=1e-12,max_step=1e-3)
nA_,TA_,UB=sol.y[:,-1]; pB=UB/1.5/2
nB_=4-nA_; TB_=pB*2/nB_
print('numeric pf',nA_*TA_,'TA',TA_,'TB',TB_,'dn',3-nA_)
assert abs(nA_*TA_-5/3)<1e-6 and abs(TB_-1.762198)<1e-5 and abs(TA_-0.790480)<1e-5; print('OK')
