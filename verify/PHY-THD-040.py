# O20223-1: 3部屋シリンダー（断熱・自由ピストン）での加熱の順序
import sympy as sp
g,V0,T0,R,Cv=sp.symbols('gamma V0 T0 R C_V',positive=True)
p0=R*T0/V0
# (1) Aのみ加熱, V_A=2V0
v=sp.symbols('v',positive=True)
sol=sp.solve(sp.Eq(2*V0+2*v,3*V0),v)[0]
p1=p0*(V0/sol)**g
TB1=p1*sol/R; TA1=p1*2*V0/R
QA1=sp.simplify(Cv*(TA1+2*TB1-3*T0))
print('(1) V_B=',sol,' T_B=',sp.simplify(TB1),' Q_A=',sp.simplify(QA1/(Cv*T0)))
# 最終状態 6:1:2
VB=V0/3; pf=p0*(V0/VB)**g
TB=sp.simplify(pf*VB/R); TA=6*TB; TC=2*TB
Q=sp.simplify(Cv*(TA+TB+TC-3*T0)); print('Q total/CvT0=',sp.simplify(Q/(Cv*T0)))
# Route I: Aを先に（第2段階でA,Bの比6:1保存）, 第1段階終わり B=C=v, A=6v
vI=sp.solve(sp.Eq(6*v+2*v,3*V0),v)[0]; pI=p0*(V0/vI)**g
UI=Cv*(pI*(6*vI+2*vI)/R)  # 全内部エネルギー
QA_I=sp.simplify(UI-3*Cv*T0); QC_I=sp.simplify(Q-QA_I)
# Route II: Cを先に（第2段階でB,Cの比1:2保存）, 第1段階終わり A=B=u, C=2u
uII=sp.solve(sp.Eq(4*v,3*V0),v)[0]; pII=p0*(V0/uII)**g
UII=Cv*(pII*4*uII/R)
QC_II=sp.simplify(UII-3*Cv*T0); QA_II=sp.simplify(Q-QC_II)
print('stage1 I: V_A=',6*vI,' II: V_C=',2*uII)
print('QA_I/CvT0=',sp.simplify(QA_I/(Cv*T0)),' QA_II/CvT0=',sp.simplify(QA_II/(Cv*T0)))
num={g:sp.Rational(5,3),Cv:sp.Rational(3,2)*R}
for name,e in [('QA_I',QA_I),('QA_II',QA_II),('QC_I',QC_I),('QC_II',QC_II),('Q',Q)]:
    print(name,'/RT0=',sp.N((e/(R*T0)).subs(num),4))
# 整合性チェック: route Iの第2段階でA,Bが断熱(pV^g一定)になっているか
VAf=2*V0  # 6/9*3V0
print('check A adiabat I:',sp.simplify((pI*(6*vI)**g - pf*VAf**g)))
print('check C adiabat II:',sp.simplify((pII*(2*uII)**g - pf*(2*VB)**g)))
print('check B adiabat:',sp.simplify(pI*vI**g-p0*V0**g), sp.simplify(pII*uII**g-p0*V0**g))
# 独立検証：数値で微小ステップの準静的加熱をシミュレート（A加熱→C加熱）
import numpy as np
from scipy.optimize import brentq
gam=5/3; Rn=1.0; T0n=1.0; V0n=1.0; Cvn=1.5
K=[1.0,1.0,1.0]  # pV^g 不変量 (p0=1, V0=1)
def state(K):
    # 共通圧力 p: sum (K_i/p)^(1/g) = 3
    f=lambda p: sum((k/p)**(1/gam) for k in K)-3
    p=brentq(f,1e-6,1e6); V=[(k/p)**(1/gam) for k in K]; return p,V
def U(K):
    p,V=state(K); return Cvn*p*sum(V)/Rn
def heat_route(order):
    K=[1.0,1.0,1.0]; Qs={0:0.0,2:0.0}
    # 目標: 最終K
    pf=3**gam; Kf=[pf*2**gam, pf*(1/3)**gam, pf*(2/3)**gam]
    for i in order:
        N=4000; k0=K[i]
        for j in range(N):
            Kold=list(K); K[i]=k0+(Kf[i]-k0)*(j+1)/N
            # dQ = dU (外部への仕事なし)
            Qs[i]+=U(K)-U(Kold)
    return Qs, state(K)
print('num route I (A then C):',heat_route([0,2]))
print('num route II (C then A):',heat_route([2,0]))
# 別法: 部屋Aだけを系にとり Q_A = ΔU_A + ∫p dV_A を数値積分
def QA_by_work(order):
    K=[1.0,1.0,1.0]; pf=3**gam; Kf=[pf*2**gam, pf*(1/3)**gam, pf*(2/3)**gam]
    QA=0.0
    for i in order:
        N=4000; k0=K[i]
        for j in range(N):
            p0_,V0_=state(K); K[i]=k0+(Kf[i]-k0)*(j+1)/N; p1_,V1_=state(K)
            dUA=Cvn*(p1_*V1_[0]-p0_*V0_[0])/Rn; W=0.5*(p0_+p1_)*(V1_[0]-V0_[0])
            QA+=dUA+W
    return QA
print('QA by work: route I',QA_by_work([0,2]),' route II',QA_by_work([2,0]))
# (4) dQ = C_V T d(pV^g)/(pV^g) の記号検証
p,V,dp,dV=sp.symbols('p V dp dV',positive=True)
dQ=(Cv/R)*(V*dp+p*dV)+p*dV
Kexpr=p*V**g
dK=sp.diff(Kexpr,p)*dp+sp.diff(Kexpr,V)*dV
print('dQ - Cv*T*dK/K =',sp.simplify(dQ-Cv*(p*V/R)*dK/Kexpr).subs(g,(Cv+R)/Cv).simplify())
