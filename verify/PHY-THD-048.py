# PHY-THD-048 排除体積をもつ気体の断熱変化
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
n,R,b,V,T,C,a,V1,V2,T1=sp.symbols('n R b V T C a V1 V2 T1',positive=True)
# dT/dV from dE=-p dV
dTdV=-(sp.Rational(2,3))*T/(V-n*b)
Tf=sp.Function('T')
sol=sp.dsolve(sp.Eq(Tf(V).diff(V),-(sp.Rational(2,3))*Tf(V)/(V-n*b)),Tf(V)); print(sol)
# power ansatz
Ta=C*(V-n*b)**a
print('a=',sp.solve(sp.simplify(sp.diff(Ta,V)-dTdV.subs(T,Ta)),a))
T2=T1*((V1-n*b)/(V2-n*b))**sp.Rational(2,3)
W=sp.Rational(3,2)*n*R*(T2-T1); print('W_on =',sp.simplify(W))
# numeric: integrate W=-∫p dV along adiabat
nn,RR,bb,v1,v2,t1=1.0,8.314,3e-5,1e-3,1e-4,300.0
f=lambda Vv,y:[-(2/3)*y[0]/(Vv-nn*bb), -nn*RR*y[0]/(Vv-nn*bb)]
s=solve_ivp(f,[v1,v2],[t1,0],rtol=1e-11,atol=1e-12)
T2n=t1*((v1-nn*bb)/(v2-nn*bb))**(2/3)
print('T2 num',s.y[0,-1],'formula',T2n,' work on gas num',-s.y[1,-1]*(-1), 'formula',1.5*nn*RR*(T2n-t1))
print('ideal T2',t1*(v1/v2)**(2/3))
# Cp check: at const p, p(V-nb)=nRT -> p dV = nR dT  => Cp=Cv+R
# (5) numbers: V1-nb = 8 (V2-nb) -> T2=4T1
print('ratio 8 ->',8**(2/3))
print('check a=-2/3:',sp.simplify(sp.diff(Ta,V)-dTdV.subs(T,Ta)).subs(a,-sp.Rational(2,3)))
# (5) after free expansion back to V1 at T=4T1, adiabatic expansion to T1: (V3-nb)=(V1-nb)*4**(3/2)
print('V3-nb =',4**1.5,'*(V1-nb)')
