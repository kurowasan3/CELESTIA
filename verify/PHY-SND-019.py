# Q125-1: 音速が高さとともに直線的に減る大気中の音線（円弧）と影の領域
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
H,h,hp,rho,y,a0=sp.symbols('H h hp rho y alpha0',positive=True)
# (1) 円の中心が高さH、半径rhoのとき cos(alpha)/c が一定であることを確認
c0=sp.Symbol('c0',positive=True)
c=c0*(1-y/H); cosal=(H-y)/rho
print('(1) cos(alpha)/c =',sp.simplify(cosal/c))
# 発射点高さ y0, 仰角 a0 -> rho=(H-y0)/cos a0
# (2) 地面に接する音線: 半径H, 中心高さH。Sを通る -> 水平距離
L1=sp.sqrt(H**2-(H-h)**2); print('(2) L=',sp.simplify(L1))
# 一般の半径rho>=Hの円が地面に達する点までの距離が rho=H で最大であることの確認
r=sp.Symbol('r',positive=True)
f=sp.sqrt(r**2-(H-h)**2)-sp.sqrt(r**2-H**2)
print('df/dr =',sp.simplify(sp.diff(f,r)))
# (3) 受け手の高さ hp
L3=sp.sqrt(h*(2*H-h))+sp.sqrt(hp*(2*H-hp))
Hv=340*1000/4.0   # 数値検証用の H（地上 340 m/s，高さ 1000 m で 336 m/s）
# --- 独立検証: 音線方程式 d/ds(n dr/ds)=grad n, n=c0/c を数値積分
Hn=Hv
def n(yy): return 1/(1-yy/Hn)
def dn(yy): return (1/Hn)/(1-yy/Hn)**2
def rhs(s,Y):
    x,yy,px,py=Y   # p = n * 単位接ベクトル
    nn=n(yy); return [px/nn,py/nn,0.0,dn(yy)]
def ground(s,Y): return Y[1]
ground.terminal=True; ground.direction=-1
def shoot(y0,alpha,smax=1500,ev=None,ms=0.5):
    nn=n(y0); return solve_ivp(rhs,[0,smax],[0,y0,nn*np.cos(alpha),nn*np.sin(alpha)],rtol=1e-10,atol=1e-10,events=ev,dense_output=True,max_step=ms)
def hit(al):
    s=shoot(2.0,al,ev=ground)
    return s.y_events[0][0][0] if len(s.t_events[0]) else None
lo,hi=-0.012,-0.004   # lo: 地面に届く, hi: 届かない
for _ in range(50):
    mid=(lo+hi)/2
    if hit(mid) is None: hi=mid
    else: lo=mid
best=hit(lo)
print('数値: 地面に届く最遠点 =',best,' 理論 L2=',float(L1.subs({H:Hn,h:2})))
assert abs(best-float(L1.subs({H:Hn,h:2})))<3
s=shoot(0.0,0.0,smax=3000,ms=5); xs=np.linspace(0,2500,25001); ys=s.sol(xs*1.0)[1]
# s≒x (ほぼ水平) なので弧長で評価してから x 座標を取る
ss=np.linspace(0,3000,300001); Ys=s.sol(ss); k=np.where(Ys[1]>=10)[0][0]
print('数値: 接点から高さ10mまで',Ys[0][k],' 理論',float(sp.sqrt(10*(2*Hn-10))))
assert abs(Ys[0][k]-float(sp.sqrt(10*(2*Hn-10))))<2
# 曲率半径の確認：(1)の円との一致
s=shoot(5.0,0.3,smax=20000,ms=20); X,Y=s.y[0],s.y[1]; rr=(Hn-5)/np.cos(0.3)
xc=0+rr*np.sin(0.3)  # 中心: 発射点から法線方向(左上)…中心 = (x0 - (-sin)…)
xc=-(-rr*np.sin(0.3))*(-1)  # 接ベクトル(cos,sin), 上側法線(-sin,cos)
xc=-rr*np.sin(0.3); yc=5+rr*np.cos(0.3)
print('中心高さ',yc,'(H=',Hn,')  円からのずれ最大',np.max(np.abs(np.hypot(X-xc,Y-yc)-rr)))
assert abs(yc-Hn)<1e-6 and np.max(np.abs(np.hypot(X-xc,Y-yc)-rr))<1e-3
# (4) 新設問：気温逆転 c=c0(1+y/K)。円の中心は高さ -K，R=(h+K)/cos(a0)
K,R_,a=sp.symbols('K R alpha',positive=True)
cinv=c0*(1+y/K); cosal2=(y+K)/R_
print('(4a) cos(alpha)/c =',sp.simplify(cosal2/cinv))
assert sp.simplify(sp.diff(cosal2/cinv,y))==0
Xf=lambda hh,KK,al: ((hh+KK)*np.sin(al)+np.sqrt((hh+KK)**2-KK**2*np.cos(al)**2))/np.cos(al)
Kn=5.0e4; hn=2.0
def n2(yy): return 1/(1+yy/Kn)
def dn2(yy): return -(1/Kn)/(1+yy/Kn)**2
def rhs2(s,Y):
    x,yy,px,py=Y; nn=n2(yy); return [px/nn,py/nn,0.0,dn2(yy)]
def ground2(s,Y): return Y[1]
ground2.terminal=True; ground2.direction=-1
prev=-1
for al in [0.0,0.002,0.005,0.01,0.02]:
    nn=n2(hn)
    sol=solve_ivp(rhs2,[0,1e5],[0,hn,nn*np.cos(al),nn*np.sin(al)],rtol=1e-11,atol=1e-11,events=ground2,max_step=20)
    Xnum=sol.y_events[0][0][0]; Xth=Xf(hn,Kn,al)
    print('(4b) alpha0=',al,' X numeric',Xnum,' formula',Xth)
    assert abs(Xnum-Xth)<1e-3*Xth and Xth>prev
    prev=Xth
assert abs(Xf(hn,Kn,0.0)-np.sqrt(hn*(hn+2*Kn)))<1e-9
# 単調増加・発散（数値）と，下向きの音線で 0〜X(0) が埋まること
als=np.linspace(0,np.pi/2-1e-6,200001); Xs=Xf(hn,Kn,als)
assert np.all(np.diff(Xs)>0) and Xs[-1]>1e9
for al in [-0.3,-0.05,-0.01]:
    nn=n2(hn)
    sol=solve_ivp(rhs2,[0,1e5],[0,hn,nn*np.cos(al),nn*np.sin(al)],rtol=1e-11,atol=1e-11,events=ground2,max_step=5)
    xd=sol.y_events[0][0][0]; print('   downward alpha0=',al,' X=',xd); assert 0<xd<Xf(hn,Kn,0.0)
print('OK')
