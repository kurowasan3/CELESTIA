# P59-1 電場の圧力と気体のアナロジー
import sympy as sp
Q,e0,a,d,x,p0,V1,nR,T1,Vb=sp.symbols('Q epsilon0 a d x p0 V1 nR T1 V_b',positive=True)
U=Q**2*d/(2*e0*a*x)
F=-sp.diff(U,x); print('(1) F=',F,' per area/(e0E^2/2):',sp.simplify(F/(a*d)/(e0*(Q/(e0*a*x))**2/2)))
pE=F/(a*d); Vol=a*d*x
print('(2) pE*V^2=',sp.simplify(pE*Vol**2),' U/(pE V)=',sp.simplify(U/(pE*Vol)))
# (3) 状態1: pE1 = p0/2, gas p = p0/2, V=V1 ; 加熱して 2V1
Vv=sp.symbols('Vv',positive=True)
pEv=p0/2*(V1/Vv)**2
pg=p0-pEv
nRT1=p0/2*V1
T2=pg.subs(Vv,2*V1)*2*V1/nR
print('(3) T2/T1=',sp.simplify(T2/(nRT1/nR)))
dU=sp.Rational(3,2)*(pg.subs(Vv,2*V1)*2*V1-nRT1)
Wg=sp.integrate(pg,(Vv,V1,2*V1))
print('   dU=',sp.simplify(dU),' Wg=',sp.simplify(Wg),' Qheat=',sp.simplify(dU+Wg))
UE1=pEv.subs(Vv,V1)*V1; UE2=pEv.subs(Vv,2*V1)*2*V1
print('   field energy decrease',sp.simplify(UE1-UE2),' = ∫pE dV',sp.simplify(sp.integrate(pEv,(Vv,V1,2*V1))))
# 別解：全体のエネルギー収支 heat = dU_gas + p0 dV - (UE1-UE2)
print('   alt heat',sp.simplify(dU+p0*V1-(UE1-UE2)))
# (4) 電池接続：pE = p0/2 一定
Qh4=sp.Rational(5,2)*p0/2*V1
# 電池の仕事: 電荷 Q ∝ x（電圧一定）→ ΔQ = Q1, W = Q1 V = 2 UE1
Wbat=2*(p0/2*V1)
dUE4=p0/2*(2*V1)-p0/2*V1
WE=p0/2*V1
print('(4) heat',Qh4,' T2/T1=',sp.Rational(1,1)*(p0/2*2*V1)/(p0/2*V1),' Wbat',Wbat,' dUE',dUE4,' work by field pressure',WE,' check',sp.simplify(Wbat-dUE4-WE))
# 電圧一定での電場の圧力を V で： pE = e0 Vb^2/(2 d^2) は x によらない
C=e0*a*x/d; Ueff=-C*Vb**2/2
print('   pE (V const)=',sp.simplify(-sp.diff(Ueff,x)/(a*d)))
