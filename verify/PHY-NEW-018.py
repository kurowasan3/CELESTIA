# PHY-NEW-018 検証：加速する台車上の粗い斜面の小物体（摩擦角と自己ロック，加速度の切り替え）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,a,th,mu=sp.symbols('g a theta mu',positive=True)
# 斜面上向き方向の力/m（静止時に必要な摩擦）: f_need = a cos - g sin ; N/m = g cos + a sin
a1=sp.solve(sp.Eq(a*sp.cos(th)-g*sp.sin(th), mu*(g*sp.cos(th)+a*sp.sin(th))),a)[0]
print('(1) a1 =',sp.simplify(a1))
phi=sp.atan(mu)
print('    tan(theta+phi) form:',sp.simplify(sp.expand_trig(g*sp.tan(th+phi))-a1))
# (2) a->oo: 限界 cos - mu sin <=0
print('(2) limit ratio', sp.limit((a*sp.cos(th)-g*sp.sin(th))/(g*sp.cos(th)+a*sp.sin(th)),a,sp.oo))
S,C=sp.Rational(28,100),sp.Rational(96,100); MU,MUK=sp.Rational(1,2),sp.Rational(1,4)
print('    numeric a1/g',(S+MU*C)/(C-MU*S), float((S+MU*C)/(C-MU*S)))
A1=2*(C-MUK*S)-(S+MUK*C); dec=S+C/2+MUK*(C-S/2); down=S+C/2-MUK*(C-S/2)
print('(3) A1',A1,' dec',dec,' s2/s',A1/dec,float(A1/dec),' static check need',S+C/2,'vs max',MU*(C-S/2),' down',down)
vP2=2*down*(1+A1/dec); print('    vP^2/(gs)=',vP2,float(vP2),' vP=',float(sp.sqrt(vP2)))
# 数値シミュレーション（g=1,s=1）：相対座標 x（斜面上向き），静止/滑りの判定込み
gn=1.0; Sn,Cn=0.28,0.96
def acc_cart(t,state):  # 台車の前向き加速度
    return 2.0 if state['phase']==1 else -0.5
def f(x,v,ac):
    N=gn*Cn+ac*Sn
    drive=ac*Cn-gn*Sn
    if abs(v)>1e-12: return drive-0.25*N*np.sign(v)
    if abs(drive)<=0.5*N: return 0.0
    return drive-0.25*N*np.sign(drive)
x,v,t,dt=0.0,0.0,0.0,1e-5; state={'phase':1}; maxx=0; vR=None
while True:
    ac=acc_cart(t,state)
    # RK2
    k1=f(x,v,ac); xm=x+0.5*dt*v; vm=v+0.5*dt*k1
    if v>0 and vm<0: vm=0
    k2=f(xm,vm,ac); vn=v+dt*k2; xn=x+dt*vm
    if v>0 and vn<0: vn=0.0
    x,v,t=xn,vn,t+dt
    if state['phase']==1 and x>=1.0: state['phase']=2; vR=v
    maxx=max(maxx,x)
    if state['phase']==2 and x<=0: break
print('sim vR^2',vR**2,' theory',float(2*A1),'  max x',maxx,' theory',float(1+A1/dec),'  vP',abs(v),' theory',float(sp.sqrt(vP2)))
