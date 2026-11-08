# Q102-1 回転する輪の張力から円弧・扇形の重心
import sympy as sp
M,R,w,a,th,s,t=sp.symbols('M R omega alpha theta s t',positive=True)
lam=M/(2*sp.pi*R)
# (1) 微小部分: T*dth = lam*R*dth*R*w^2
T=lam*R**2*w**2; print('T',sp.simplify(T))
# (2) 円弧の重心（直接積分）
dc=sp.integrate(R*sp.cos(th),(th,-a,a))/(2*a); print('arc COM',sp.simplify(dc))
# 張力の合力 2T sin a = (2a/2pi)M d w^2 から
d2=sp.solve(sp.Eq(2*T*sp.sin(a),(2*a/(2*sp.pi))*M*sp.Symbol('d')*w**2),sp.Symbol('d'))[0]
print('from tension',sp.simplify(d2-dc))
# (3) 切断後
v1=w*R*sp.sin(a)/a; v2=w*R*sp.sin(sp.pi-a)/(sp.pi-a)
m1=M*a/sp.pi; m2=M*(sp.pi-a)/sp.pi
print('momentum',sp.simplify(m1*v1-m2*v2))
print('ratio v1/v2',sp.simplify(v1/v2))
D0=R*sp.sin(a)*(1/a+1/(sp.pi-a)); Dt=sp.sqrt(D0**2+((v1+v2)*t)**2)
print('D(t)/D0',sp.simplify(Dt/D0))
# 数値: 実際に各断片の点を自由運動させて重心距離を計算
import numpy as np
Rv,wv,av,tv=1.0,2.0,0.7,1.3
def piece(lo,hi,n=4000):
    th_=np.linspace(lo,hi,n)
    x=Rv*np.cos(th_); y=Rv*np.sin(th_)
    vx=-wv*y; vy=wv*x
    return (x+vx*tv).mean(),(y+vy*tv).mean()   # 各点が等速直線運動するわけではないが重心は等速
# 重心は等速直線運動なので初期重心+初期重心速度*t で計算
def com(lo,hi,n=200001):
    th_=np.linspace(lo,hi,n); x=Rv*np.cos(th_);y=Rv*np.sin(th_)
    return np.array([x.mean(),y.mean()]),np.array([(-wv*y).mean(),(wv*x).mean()])
c1,u1=com(-av,av); c2,u2=com(av,2*np.pi-av)
P1=c1+u1*tv; P2=c2+u2*tv
print('D num',np.linalg.norm(P1-P2),' theory',float(Dt.subs({R:Rv,w:wv,a:av,t:tv})))
# (4) 扇形の重心
ds=sp.integrate(sp.integrate(s*sp.cos(th)*s,(s,0,R)),(th,-a,a))/sp.integrate(sp.integrate(s,(s,0,R)),(th,-a,a))
print('sector',sp.simplify(ds), ' half disk', sp.simplify(ds.subs(a,sp.pi/2)))
