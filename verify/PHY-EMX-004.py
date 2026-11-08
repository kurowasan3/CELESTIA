# Q187-1 磁場中を流れる誘電体流体と極板の電圧・電荷
import sympy as sp
e0,er,v,B,d,S,E,q,w,l,t,C=sp.symbols('epsilon0 epsilon_r v B d S E q w ell t C',positive=True)
# 方法1：分極電荷の面密度モデル（E は上→下向きを正, E'=vB-E 上向き）
sig=e0*(er-1)*(v*B-E)
# (1) 開放：極板電荷0, 極板間の電場は分極電荷だけでつくられる
E1=sp.solve(sp.Eq(E,sig/e0),E)[0]; V1=sp.simplify(E1*d); print('V_open=',sp.factor(V1))
# 方法2：電束密度 D = e0 E_lab + P, 自由電荷0なら D=0（符号：上向き正で）
Eu=sp.symbols('Eu'); P=e0*(er-1)*(Eu+v*B)   # 上向き成分
Eu_s=sp.solve(sp.Eq(e0*Eu+P,0),Eu)[0]; assert sp.simplify(-Eu_s*d-V1)==0
# (2) 短絡：V=0 → E=0, 上の極板の電荷 = -sig*S
Q2=sp.simplify(-sig.subs(E,0)*S); print('Q_upper(short)=',Q2)
# (3) 流入：入った部分の面積 w*v*t → 電荷の大きさ e0(er-1)vB w v t, 電流 = d/dt
I3=sp.diff(e0*(er-1)*v*B*w*v*t,t); print('I_fill=',I3)
# 外部コンデンサーC：上極板電荷 q, E=(q/S+sig)/e0, 外部Cの電圧 = -q/C
Eq=sp.solve(sp.Eq(E,(q/S+sig)/e0),E)[0]; qq=sp.solve(sp.Eq(Eq*d,-q/C),q)[0]
Vc=sp.simplify(-qq/C); C0=er*e0*S/d
assert sp.simplify(Vc - V1*C0/(C0+C))==0; print('V with C =',sp.factor(Vc))
print('limit er->oo:',sp.limit(V1,er,sp.oo))
print('OK')
