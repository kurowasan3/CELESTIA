# P26 類題：浮かぶ半球殻の縁におもりをのせる（浮力の作用線は球の中心を通る）
import numpy as np
from scipy.optimize import brentq
R=1.0; rho=1.0; g=1.0
m=5*np.pi*rho*R**3/24
Vcap=lambda d: np.pi*d**2*(3*R-d)/3
assert abs(rho*Vcap(R/2)-m)<1e-12
def state(k):
    th=np.arctan(2*k)
    d=brentq(lambda d: rho*Vcap(d)-m*(1+k),0,R)
    hA=(R-d)-R*np.sin(th)
    return th,d,hA
for k in [0.2]:
    th,d,hA=state(k); print('k=',k,' theta(deg)=',np.degrees(th),' d/R=',d,' hA/R=',hA)
kc=brentq(lambda k: state(k)[2],0.01,1.0); print('critical m\'/m =',kc, ' theta=',np.degrees(np.arctan(2*kc)),' d=',state(kc)[1])
# 独立検証：傾いた半球殻の表面に水圧を直接積分して，力とOまわりのトルクを計算
def pressure_force(th,d,n=800):
    # 殻の軸（Oから底へ向かう単位ベクトル）を xz 面内で鉛直下向きから +x 側へ θ 傾ける
    ax=np.array([np.sin(th),0,-np.cos(th)])
    e1=np.array([np.cos(th),0,np.sin(th)]); e2=np.array([0,1,0])
    zw=-(R-d)   # 水面の高さ（O を原点）
    F=np.zeros(3); Tq=np.zeros(3)
    a=np.linspace(0,np.pi/2,n+1); b=np.linspace(0,2*np.pi,2*n+1)
    am=(a[1:]+a[:-1])/2; bm=(b[1:]+b[:-1])/2; da=a[1]-a[0]; db=b[1]-b[0]
    A,B=np.meshgrid(am,bm,indexing='ij')
    nx=(np.cos(A)[...,None]*ax+np.sin(A)[...,None]*(np.cos(B)[...,None]*e1+np.sin(B)[...,None]*e2))
    P=R*nx; z=P[...,2]
    p=np.where(z<zw,rho*g*(zw-z),0.0)
    dS=R**2*np.sin(A)*da*db
    dF=-(p*dS)[...,None]*nx   # 外向き法線 nx，圧力は内向き
    F=dF.sum(axis=(0,1)); Tq=np.cross(P,dF).sum(axis=(0,1))
    return F,Tq
k=0.2; th,d,hA=state(k)
F,Tq=pressure_force(th,d)
print('buoyancy (numeric) Fz=',F[2],' expected',m*(1+k)*g,'  torque about O:',Tq)
# 重力のトルク：殻の重心は軸上 R/2，おもりは縁 A（-x 側で最も低い縁）
G=0.5*R*np.array([np.sin(th),0,-np.cos(th)])
Arim=R*np.array([-np.cos(th),0,-np.sin(th)])
tq=np.cross(G,[0,0,-m*g])+np.cross(Arim,[0,0,-k*m*g])
print('gravity torque y:',tq[1])
assert abs(F[2]-m*(1+k))<2e-3 and abs(Tq[1])<2e-3 and abs(tq[1])<1e-12
