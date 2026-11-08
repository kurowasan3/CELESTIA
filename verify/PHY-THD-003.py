# M51-1: sudden release of piston (irreversible adiabatic) vs quasi-static
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
P0,S,L,M,g,T0,n,R,h,T=sp.symbols('P0 S L M g T0 n R h T',positive=True)
P1=P0+M*g/S
eqs=[sp.Eq(P1*S*h,n*R*T),sp.Eq(sp.Rational(3,2)*n*R*(T-T0),(P0*S+M*g)*(L-h)),sp.Eq(P0*S*L,n*R*T0)]
sol=sp.solve(eqs,[h,T,n],dict=True)[0]
hs=sp.simplify(sol[h]); Ts=sp.simplify(sol[T]); print('h',hs,' T',Ts)
assert sp.simplify(hs-L*(5*P0*S+2*M*g)/(5*(P0*S+M*g)))==0
assert sp.simplify(Ts-T0*(1+2*M*g/(5*P0*S)))==0
a=sp.symbols('alpha',positive=True)
print('limit h/L',sp.limit((hs/L).subs(M,a*P0*S/g),a,sp.oo))
# quasi-static
hq=L*(P0/P1)**sp.Rational(3,5); Tq=T0*(P1/P0)**sp.Rational(2,5)
# return heat
Q=sp.simplify(sp.Rational(5,2)*P1*S*(L-hs)); print('Q return',Q)
# compare T >= Tq numerically
for al in [0.1,1,5,50]:
    print(al, 'T/T0',1+0.4*al,'Tq/T0',(1+al)**0.4,'h/L',(5+2*al)/(5+5*al),'hq/L',(1+al)**-0.6)
# dynamic simulation: piston with damping proportional to velocity (dissipated into gas)
# units: P0=1,S=1,L=1,nRT0=1, M g = alpha; gas energy U=1.5 P V; dissipation heats gas
al=2.0; Mg=al; Mass=al/9.8; c=0.5
def rhs(t,y):
    x,v,U=y   # x height, U internal energy
    P=U/(1.5*x)
    F=P-1.0-Mg-c*v
    dU=-P*v + c*v*v   # adiabatic work + dissipated heat
    return [v,F/Mass,dU]
sol=solve_ivp(rhs,[0,200],[1.0,0.0,1.5],rtol=1e-10,atol=1e-12,max_step=0.05)
x,v,U=sol.y[:,-1]
print('sim h',x,'formula',(5+2*al)/(5+5*al),' T/T0 sim',U/1.5,'formula',1+0.4*al,' min x',sol.y[0].min())
