# J003-1 検証：粗い斜面を上って戻る物体の通過時刻から加速度・初速・傾き・摩擦係数を求める逆問題
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
d,t1,t2,t3=sp.Rational(3,4),sp.Rational(1,4),sp.Integer(1),sp.Rational(3,2)
k=t1/(t3-t2); print('k=sqrt(a_d/a_u)=',k)
tau3=t3/(1+k); tau2=t2-k*tau3; print('tau3,tau2',tau3,tau2)
ad=2*d/(tau3**2-tau2**2); au=ad/k**2; print('a_d, a_u =',ad,au)
v0=au*k*tau3; print('v0 =',v0)
g=sp.Rational(98,10)
s=(au+ad)/(2*g); muc=(au-ad)/(2*g); c=sp.sqrt(1-s**2); mu=muc/c
print('sin th =',s,sp.N(s),' theta deg',sp.N(sp.asin(s)*180/sp.pi),' mu =',sp.N(mu))
# 数値積分で通過時刻を再現
au_,ad_=float(au),float(ad); v0_=float(v0)
def f(t,y):
    x,v=y
    return [v, -au_ if v>0 else ad_*(-1)]
# 上り: 加速度 -au, 下り: 加速度 -ad（x 正が上向き）
def rhs(t,y):
    x,v=y; a=-au_ if v>0 else -ad_
    return [v,a]
ev=[]
for target in [0.75,0.0]:
    pass
sol=solve_ivp(rhs,[0,2.0],[0,v0_],max_step=1e-4,dense_output=True)
ts=np.linspace(0,1.6,1600001); xs=sol.sol(ts)[0]
cross=ts[np.where(np.diff(np.sign(xs-0.75))!=0)[0]]; print('times at x=0.75:',cross)
cross0=ts[1:][np.where(np.diff(np.sign(xs[1:]-0.0))!=0)[0]]; print('times at x=0:',cross0)
# 同じ点での上りと下りの速さの比 = sqrt(a_u/a_d)
H=v0_**2/(2*au_); x=0.3
print('speed ratio at x=0.3:',np.sqrt(v0_**2-2*au_*x)/np.sqrt(2*ad_*(H-x)),' sqrt(au/ad)=',np.sqrt(au_/ad_))
