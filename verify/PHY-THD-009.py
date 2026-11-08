import sympy as sp, numpy as np
n,R,T0,V0,VB,u=sp.symbols('n R T0 V0 V_B u',positive=True)
P0=n*R*T0/V0
# energy: U_B = (3/2) n' R T0 + P0 u, n' = P0 u/(R T0)
UB=sp.Rational(3,2)*(P0*u/(R*T0))*R*T0+P0*u
PB=sp.simplify(UB/(sp.Rational(3,2)*VB)); print('PB',PB, ' ratio PB/P0 =',sp.simplify(PB/P0))
print('VB min',sp.solve(sp.Eq(PB.subs(u,V0),P0),VB))
Tf=sp.simplify((UB.subs(u,V0))/(sp.Rational(3,2)*n*R)); print('Tf',Tf)
# (3) B initially n, T0
UB3=sp.Rational(3,2)*n*R*T0+sp.Rational(5,2)*P0*u
Tf3=sp.simplify(UB3.subs(u,V0)/(sp.Rational(3,2)*2*n*R)); print('Tf3',Tf3)
PB3=sp.simplify(UB3/(sp.Rational(3,2)*VB)); print('PB3',PB3)
print('VB min3',sp.solve(sp.Eq(PB3.subs(u,V0),P0),VB))
# independent numeric: discrete parcels; each parcel dn at T0 pushed in with work P0*dV, B mixes
N=200000; Vb=2.0; nn=1.0; RR=1; TT=1; VV=1; P=nn*RR*TT/VV
U=0.0; nb=0.0; du=VV/N
for i in range(N):
    dn=P*du/(RR*TT); U+=1.5*dn*RR*TT+P*du; nb+=dn
print('numeric Tf',U/(1.5*nb*RR),'PB/P0',U/1.5/Vb/P)
assert abs(U/(1.5*nb)-5/3)<1e-9
U=1.5; nb=1.0; Vb=3.0
for i in range(N):
    dn=P*du; U+=1.5*dn+P*du; nb+=dn
print('numeric Tf3',U/(1.5*nb)); assert abs(U/(1.5*nb)-4/3)<1e-9; print('OK')
