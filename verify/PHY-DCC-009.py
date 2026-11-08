# P63-1 テーパーのある抵抗線のメートルブリッジと入れ替え法
import sympy as sp, math
from scipy.integrate import quad
rho,a,L,x,beta,Rs,Rx=sp.symbols('rho a L x beta R_s R_x',positive=True)
s=sp.symbols('s',positive=True)
r=a*(1+(beta-1)*s/L)
RAD=sp.simplify(sp.integrate(rho/(sp.pi*r**2),(s,0,x)))
RAB=sp.simplify(RAD.subs(x,L))
print('(1) R_AD=',RAD,'  R_AB=',RAB)
ratio=sp.simplify(RAD/(RAB-RAD)); print('   R_AD/R_DB=',ratio)
# (2) 数値: 1回目 Rs(AC)/Rx(CB) = ratio(x1); 入れ替え Rx/Rs = ratio(x2)
x1,x2=0.250,0.571
b2=(1-x1)*(1-x2)/(x1*x2); b=math.sqrt(b2)
Rx_=10.0*(1-x1)/(b*x1)
print('(2) beta=',b,' Rx=',Rx_, ' alt from x2:',10.0*b*x2/(1-x2))
# 真値 beta=1.5, Rx=20 から読みを再現
bt,Rxt=1.5,20.0
x1t=sp.nsolve(bt*x/(1-x)-10/Rxt,x,0.3); x2t=sp.nsolve(bt*x/(1-x)-Rxt/10,x,0.5)
print('   readings from truth',x1t,x2t)
# 数値積分で R_AD 比を確認
f=lambda s_: 1/(1+(bt-1)*s_)**2
print('   numeric ratio at x1',quad(f,0,0.25)[0]/quad(f,0.25,1)[0],' expected',bt*0.25/0.75)
# (3) 一様と仮定した場合
naive1=10*(1-x1)/x1; naive2=10*x2/(1-x2)
print('(3) naive1',naive1,' naive swap geometric mean',math.sqrt(naive1*naive2))
