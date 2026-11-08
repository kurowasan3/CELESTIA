# Soap bubble (2 surfaces, surface tension sigma) containing gas A inside cylinder of gas B at const P
import sympy as sp
P,s,r,R,nA,nB,T=sp.symbols('P sigma r R n_A n_B T',positive=True)
PA=P+4*s/r
VA=sp.Rational(4,3)*sp.pi*r**3
# force balance on hemisphere: (PA-P)*pi r^2 = 2 surfaces * sigma * 2 pi r
print('excess p:',sp.solve(sp.Eq((sp.Symbol('dp'))*sp.pi*r**2,2*s*2*sp.pi*r),sp.Symbol('dp')))
# membrane energy: (PA-P) dVA vs d(sigma*2*4pi r^2)
print('film energy check:',sp.simplify((PA-P)*sp.diff(VA,r)-sp.diff(s*8*sp.pi*r**2,r)))
# T as function of r (gas A eq of state)
TA=PA*VA/(nA*R)
dTdr=sp.diff(TA,r)
# heat: dQ = 3/2(nA+nB) R dT + P dVB + PA dVA ; P VB = nB R T
dQdr=sp.Rational(3,2)*(nA+nB)*R*dTdr+nB*R*dTdr+PA*sp.diff(VA,r)
C=sp.simplify(dQdr/dTdr)
target=sp.Rational(5,2)*nB*R+nA*R*(sp.Rational(3,2)+(3*P*r+12*s)/(3*P*r+8*s))
print('C matches:',sp.simplify(C-target)==0)
print('limit sigma->0:',sp.simplify(sp.limit(C,s,0)),' limit P->0:',sp.simplify(sp.limit(C,P,0)))
# independent check: energy bookkeeping via total enthalpy-like route: Q = dU + P dV_total + d(film energy)
VB=nB*R*T/P
# param by r: T(r)
Tr=TA
Utot=sp.Rational(3,2)*(nA+nB)*R*Tr
Vtot=VA+VB.subs(T,Tr)
Efilm=s*8*sp.pi*r**2
dQ2=sp.diff(Utot,r)+P*sp.diff(Vtot,r)+sp.diff(Efilm,r)
print('route2 equal:',sp.simplify(dQ2-dQdr)==0)
# 参考（数値例）: r=1.0e-2 m, sigma=2.5e-2, P=1.0e5: relative deviation of A's molar heat capacity from 5/2 R
x=(3*P*r+12*s)/(3*P*r+8*s)
val=((sp.Rational(3,2)+x)/sp.Rational(5,2)-1).subs({P:1.0e5,r:1.0e-2,s:2.5e-2})
print('relative deviation:',sp.N(val),' approx 8 sigma/(15 P r)=',sp.N((8*s/(15*P*r)).subs({P:1.0e5,r:1.0e-2,s:2.5e-2})))
# 参考: radius doubling: T2/T1 with x=4 sigma/(P r1)
r1,xx=sp.symbols('r1 x',positive=True)
ratio=sp.simplify(TA.subs(r,2*r1)/TA.subs(r,r1))
print('T2/T1=',sp.simplify(ratio.subs(s,xx*P*r1/4)))
# heat for doubling radius expressed
# ---- 新しい (5)：外が真空（P=0）のとき P_A V_A^k = 一定，k=1/3，ポリトロープのモル比熱
PA0=4*s/r
kk=sp.simplify(-sp.diff(sp.log(PA0),r)/sp.diff(sp.log(VA),r)); print('(5a) k=',kk); assert kk==sp.Rational(1,3)
assert sp.simplify(sp.diff(PA0*VA**sp.Rational(1,3),r))==0
k_=sp.symbols('k',positive=True); V_=sp.symbols('V',positive=True); n_=sp.symbols('n',positive=True)
Pp=1/V_**k_                        # P V^k = 1
dW=Pp; dnRT=sp.diff(Pp*V_,V_)      # 単位体積変化あたり
Cpoly=sp.simplify((sp.Rational(3,2)*dnRT+dW)/dnRT*R)
print('(5b) C=',Cpoly); assert sp.simplify(Cpoly-(sp.Rational(3,2)*R+R/(1-k_)))==0
assert sp.simplify(Cpoly.subs(k_,sp.Rational(1,3))-3*R)==0 and sp.simplify(Cpoly.subs(k_,0)-sp.Rational(5,2)*R)==0
cA=sp.simplify((C-sp.Rational(5,2)*nB*R)/nA)
assert sp.simplify(cA.subs(P,0)-3*R)==0
print('    k=1/3 → 3R，(3) の式で P=0 としても c_A=3R')
print('OK (5)')
