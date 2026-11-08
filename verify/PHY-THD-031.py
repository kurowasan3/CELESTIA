# J068-1 検証：クレマン・デゾルムの方法（断熱膨張＋定積での温度回復）
import sympy as sp
p0,p1,T0,g_,V,k=sp.symbols('p0 p1 T0 gamma V k',positive=True)
# 残留気体：初め (p1, Vr, T0) -> 断熱 -> (p0, V, T') -> 定積 -> (p2, V, T0)
Vr=V*(p0/p1)**(1/g_)
Tp=T0*(p0/p1)**((g_-1)/g_)
# 状態方程式で T' を確認: p0 V / T' = p1 Vr / T0
print('eqn of state check:',sp.simplify(p0*V/Tp - p1*Vr/T0))
p2=p0*T0/Tp
print('p2 =',sp.simplify(p2))
# (2) 残留物質量比 = Vr/V (等温・同圧p1での体積比) と p2/p1 の一致
print('n ratio - p2/p1 =',sp.simplify(sp.powsimp(Vr/V - p2/p1,force=True)))
# (3) gamma = ln(p1/p0)/ln(p1/p2)
expr=sp.log(p1/p0)/sp.log(p1/p2)
print('gamma check:',sp.simplify(sp.expand_log(expr,force=True)))
# 近似: p1=p0+a h1, p2=p0+a h2 ; 1次近似
a,h1,h2,eps=sp.symbols('a h1 h2 epsilon',positive=True)
ap=sp.log((p0+eps*a*h1)/p0)/sp.log((p0+eps*a*h1)/(p0+eps*a*h2))
print('approx limit:',sp.limit(ap,eps,0))
print('numeric gamma approx:',15.0/(15.0-4.3))
# 正確な式との比較（水柱 ρ=1000, g=9.8, p0=1.013e5）
P0=1.013e5; P1=P0+1000*9.8*0.150; P2=P0+1000*9.8*0.043
import math
print('exact gamma from numbers:',math.log(P1/P0)/math.log(P1/P2))
# (4) ポリトロープ pV^k: p2 = p0 (p1/p0)^((k-1)/k) -> 式は k を与える
p2k=p0*(p1/p0)**((k-1)/k)
print('poly estimate:',sp.simplify(sp.expand_log(sp.log(p1/p0)/sp.log(p1/p2k),force=True)))
print('k->1: p2 =',sp.limit(p2k,k,1))
# (5) Q = V(p2-p0)/(gamma-1)
gam=15.0/10.7; Q=1e-2*1000*9.8*0.043/(gam-1); print('Q =',Q,' (gamma=1.40:',1e-2*1000*9.8*0.043/0.40,')')
# 別解：n Cv dT  with n = p2 V/(R T0), Cv=R/(gamma-1), dT = T0 - T' = T0(1 - p0/p2)
R=sp.symbols('R',positive=True); p2s=sp.symbols('p2',positive=True)
Qalt=(p2s*V/(R*T0))*(R/(g_-1))*(T0-T0*p0/p2s)
print('Q alt =',sp.simplify(Qalt))
