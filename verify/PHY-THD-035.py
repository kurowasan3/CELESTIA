# O20093-1 : 重いピストンで仕切られた鉛直断熱シリンダー内の熱の移動
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
R,T0,V0=1,1,1           # 単位系 R=T0=V0=1
PM=3                    # Mg/S = 3 R T0/V0
TB=1; VA0=1; VB0=3; Vt=4
pB=sp.Rational(TB,VB0); pA=pB+PM; TA=pA*VA0; print('(1) TA=',TA)
u,T=sp.symbols('u T',positive=True)
E0=sp.Rational(3,2)*(TA+TB)+PM*VA0
sol=sp.solve([sp.Eq(T*(1/u-1/(Vt-u)),PM),sp.Eq(3*T+PM*u,E0)],[u,T],dict=True)
print('(2) solutions',sol)
print('average',(TA+TB)/2,' piston drop (S*dx)=',VA0-sp.Rational(2,3),' PE released=',PM*(VA0-sp.Rational(2,3)))
# (3) ピストン固定で温度が等しくなった後: T'=(TA+TB)/2, 圧力差
Tp=(TA+TB)/2; print('fixed: T=',Tp,' pA-pB=',Tp/VA0-Tp/VB0,' vs PM=',PM)
# 手で支えながらゆっくり: 両気体共通温度, 3 dT = -(pA-pB) dVA = -T(1/VA-1/(Vt-VA)) dVA
def f(VA,y): T_=y[0]; return [-T_*(1/VA-1/(Vt-VA))/3]
def ev(VA,y): return y[0]*(1/VA-1/(Vt-VA))-PM
ev.terminal=True
s=solve_ivp(f,(1.0,0.01),[float(Tp)],events=ev,rtol=1e-12,atol=1e-12)
print('hand-lowered final VA=',s.t_events[0][0],' T=',s.y_events[0][0][0],' (< 5/2)')
# 突然放した場合: 減衰振動の数値モデル（ピストン速度に比例する内部散逸→気体へ熱として戻る, 温度は熱伝導で緩和）
m_p=1.0; S=1.0; g=PM*S/m_p
def rhs(t,y):
    x,vx,TA_,TB_=y         # x = VA/S
    VA=S*x; VB=Vt-VA
    pa=TA_/VA; pb=TB_/VB
    damp=2.0*vx            # 散逸力
    F=(pa-pb)*S-m_p*g-damp
    kappa=0.05             # ピストンを通るゆっくりした熱伝導
    q=kappa*(TA_-TB_)
    # 第一法則: 1.5 dT_A = -pa dVA - q + 散逸の半分 ; B も同様
    dTA=(-pa*S*vx - q + 0.5*damp*vx)/1.5
    dTB=( pb*S*vx + q + 0.5*damp*vx)/1.5
    return [vx,F/m_p,dTA,dTB]
y0=[1.0,0.0,float(Tp),float(Tp)]
s2=solve_ivp(rhs,(0,600),y0,rtol=1e-9,atol=1e-11,max_step=0.5)
print('sudden release final VA=',s2.y[0,-1],' TA=',s2.y[2,-1],' TB=',s2.y[3,-1])
y0=[1.0,0.0,float(TA),float(TB)]
s3=solve_ivp(rhs,(0,600),y0,rtol=1e-9,atol=1e-11,max_step=0.5)
print('free from start final VA=',s3.y[0,-1],' TA=',s3.y[2,-1],' TB=',s3.y[3,-1])
