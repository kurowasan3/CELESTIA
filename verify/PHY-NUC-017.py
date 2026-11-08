# J148-2 検証: 2種のウラン同位体の存在比・放射能比の時間変化
import sympy as sp, numpy as np
T5,T8=0.704,4.47   # 10^9 年
r0=0.0073
Teff=T5*T8/(T8-T5); print('ratio half-life (Gyr)=',Teff)
D=1/T5-1/T8; print('Delta=',D)
t_oklo=np.log2(0.030/r0)/D; print('t(3.0%)=',t_oklo)
A_now=r0*T8/T5; print('activity ratio now=',A_now)
t_eq=np.log2(T5/(T8*r0))/D; print('t(A equal)=',t_eq)
# 常用対数での手計算
l2,l3,l73,l704,l447=0.301,0.477,0.863,0.848,0.650
print('hand t_oklo=',((l3-2)-(l73-3))/l2/D, ' hand t_eq=',((l704+8)-(l447+9)-(l73-3))/l2/D)
# 直接シミュレーション: 現在 N8=1, N5=r0 から過去へ
t=np.linspace(0,4,400001)
N5=r0*2**(t/T5); N8=2**(t/T8)
print('sim t(3%)=',t[np.argmax(N5/N8>=0.030)],' sim t(Aeq)=',t[np.argmax(N5/T5>=N8/T8)])
t_s=sp.symbols('t'); print(sp.simplify(sp.Symbol('r0')*2**(t_s/sp.Symbol('T5'))/2**(t_s/sp.Symbol('T8'))))
