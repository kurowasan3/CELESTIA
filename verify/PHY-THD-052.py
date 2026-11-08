# PHY-THD-052 検証: タンク付きシリンダーの熱機関（等圧での混合は効率を下げない）
import sympy as sp
b,k=sp.symbols('beta k',positive=True)
P0,V0,T0,R=sp.symbols('P0 V0 T0 R',positive=True)
n=P0*V0/(R*T0)
s=b**sp.Rational(3,5)
# A->B 定積加熱
Qin=sp.Rational(3,2)*n*R*(b*T0-T0)
# B->C 断熱膨張 P: bP0->P0
VC=V0*b**sp.Rational(3,5); TC=T0*VC/V0
WBC=(b*P0*V0-P0*VC)/sp.Rational(2,3)
# C->D 混合（体積固定、断熱）: 内部エネルギー保存
TD=sp.symbols('TD')
TDs=sp.solve(sp.Eq(sp.Rational(3,2)*R*(n*TC+k*n*T0),sp.Rational(3,2)*R*(1+k)*n*TD),TD)[0]
PD=(1+k)*n*R*TDs/(VC+k*V0)
print('TD/T0=',sp.simplify(TDs/T0),' PD/P0=',sp.simplify(PD/P0))
ncylD=P0*VC/(R*TDs)
print('moved tank->cyl /n =',sp.simplify((ncylD-n)/n))
# D->A 定圧冷却（全体）
Qout=sp.Rational(5,2)*(1+k)*n*R*(TDs-T0)
WDA=P0*((1+k)*V0-(VC+k*V0))
W=sp.simplify(WBC+WDA)
print('Qout=',sp.simplify(Qout),'  Qin-Qout-W=',sp.simplify(Qin-Qout-W))
eta=sp.simplify(1-Qout/Qin); print('eta=',eta,' eta(b=32)=',eta.subs(b,32))
# 変形: B でコックを開く（圧力が異なる混合）→ 全体の断熱膨張 → 定圧冷却
Pm=(b+k)/(1+k)*P0
Tm=sp.solve(sp.Eq(b*T0+k*T0,(1+k)*sp.Symbol('Tm')),sp.Symbol('Tm'))[0]
print('Pm check (PV=nRT):',sp.simplify(Pm*(1+k)*V0-(1+k)*n*R*Tm))
Vf=(1+k)*V0*(Pm/P0)**sp.Rational(3,5)
Qout2=sp.Rational(5,2)*P0*(Vf-(1+k)*V0)
eta2=1-Qout2/Qin
f=sp.lambdify((b,k),eta2); g=sp.lambdify((b,),eta)
for bb in [32,243,10]:
    print(bb,[ (kk,round(f(bb,kk),5)) for kk in [1e-9,0.5,1,3,100,1e6]],' eta=',g(bb))
print('limit k->oo:',sp.limit(eta2.subs(b,32),k,sp.oo),' k->0:',sp.simplify(sp.limit(eta2,k,0)-eta))
# 別解: エンタルピー H=5/2 PV の保存（等圧混合）で Qout 確認
print('Qout via enthalpy:',sp.simplify(sp.Rational(5,2)*P0*(VC+k*V0)-sp.Rational(5,2)*P0*(1+k)*V0 - Qout))
