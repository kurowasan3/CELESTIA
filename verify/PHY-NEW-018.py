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
# --- (3) 地面から見た動摩擦力の仕事（a>a1 一定，斜面に沿って s のぼるまで）
mup,m_,s_=sp.symbols("mu' m s",positive=True)
N=m_*(g*sp.cos(th)+a*sp.sin(th)); f=mup*N
A=a*(sp.cos(th)-mup*sp.sin(th))-g*(sp.sin(th)+mup*sp.cos(th))
t2=2*s_/A; X=a*t2/2
r=sp.Matrix([X-s_*sp.cos(th), s_*sp.sin(th)])          # 地面から見た小物体の変位（前方 x, 上 y）
fvec=f*sp.Matrix([sp.cos(th),-sp.sin(th)])               # 斜面下向き
W=sp.simplify(fvec.dot(r)); Wp=sp.simplify(-fvec.dot(sp.Matrix([X,0])))
W_claim=mup*N*(g*(sp.sin(th)+mup*sp.cos(th))+mup*a*sp.sin(th))/A*s_
Wp_claim=-mup*N*a*sp.cos(th)/A*s_
print('W check',sp.simplify(W-W_claim),' Wp check',sp.simplify(Wp-Wp_claim),' sum',sp.simplify(W+Wp))
assert sp.simplify(W-W_claim)==0 and sp.simplify(Wp-Wp_claim)==0 and sp.simplify(W+Wp+mup*N*s_)==0
# 独立な確認：地面系の仕事とエネルギーの関係 W_f = ΔK + mg s sinθ - W_N（W_N は垂直抗力の仕事）
t=sp.sqrt(t2); vg=sp.Matrix([a*t-A*t*sp.cos(th), A*t*sp.sin(th)])
dK=m_*vg.dot(vg)/2; nvec=N*sp.Matrix([sp.sin(th),sp.cos(th)]); WN=nvec.dot(r)
print('energy-theorem check',sp.simplify(dK+m_*g*s_*sp.sin(th)-WN-W_claim))
assert sp.simplify(dK+m_*g*s_*sp.sin(th)-WN-W_claim)==0
# 数値：ランダムなパラメータで A>0 のとき W>0 を確認，台車系で数値積分した小物体の運動から仕事を直接積分
rng=np.random.default_rng(1)
for _ in range(2000):
    thv=rng.uniform(0.05,1.3); muv=rng.uniform(0.05,1.0); mupv=muv*rng.uniform(0.05,0.99)
    if muv*np.tan(thv)>=1: continue
    a1v=(np.sin(thv)+muv*np.cos(thv))/(np.cos(thv)-muv*np.sin(thv))
    av=a1v*rng.uniform(1.001,5)
    vals={g:1,th:thv,a:av,mup:mupv,m_:1,s_:1}
    Av=float(A.subs(vals)); assert Av>0
    assert float(W_claim.subs(vals))>0
dt=1e-5; thv,av,mupv=0.4,3.0,0.3; Nv=np.cos(thv)+av*np.sin(thv); Av=av*(np.cos(thv)-mupv*np.sin(thv))-(np.sin(thv)+mupv*np.cos(thv))
xr=0.0; u=0.0; Vc=0.0; Wnum=0.0
while xr<1.0:
    vb=np.array([Vc-u*np.cos(thv), u*np.sin(thv)])          # 地面から見た小物体の速度
    Wnum+=mupv*Nv*np.dot([np.cos(thv),-np.sin(thv)],vb)*dt
    u+=Av*dt; xr+=u*dt-Av*dt*dt/2; Vc+=av*dt
print('numeric W',Wnum,' formula',float(W_claim.subs({g:1,th:thv,a:av,mup:mupv,m_:1,s_:1})))
assert abs(Wnum-float(W_claim.subs({g:1,th:thv,a:av,mup:mupv,m_:1,s_:1})))<1e-3
print('(3) OK')
