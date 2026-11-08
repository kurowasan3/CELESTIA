# P1 類題：直線上を等速で進むBを追うA（u r + v s の時間変化が一定）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
u,v,l,th,t=sp.symbols('u v ell theta t',positive=True)
drdt=-u+v*sp.cos(th); dsdt=v-u*sp.cos(th)
print('d(ur+vs)/dt =',sp.simplify(u*drdt+v*dsdt))
T=u*l/(u**2-v**2); print('T=',T,' dist=',sp.simplify(u*T))
# 一般の初期角 alpha (Bの速度とB→Aの逆…) s0 = -l cos? ここでは s0=0
# u=v: r+s=l, 終局 theta->0 → r=s → r=l/2
# 数値
def run(uu,vv,ll=1.0,tmax=50):
    def f(t,y):
        xa,ya=y; xb,yb=vv*t,0.0
        dx,dy=xb-xa,yb-ya; r=np.hypot(dx,dy)
        return [uu*dx/r,uu*dy/r]
    ev=lambda t,y: np.hypot(vv*t-y[0],y[1])-1e-6
    ev.terminal=True
    return solve_ivp(f,[0,tmax],[0,-ll],events=ev,rtol=1e-11,atol=1e-13,max_step=1e-3)
for uu,vv in [(2.0,1.0),(3.0,2.0),(1.5,1.2)]:
    s=run(uu,vv)
    Tn=s.t_events[0][0] if len(s.t_events[0]) else None
    print(uu,vv,'T num',Tn,' theory',uu/(uu**2-vv**2),' catch x',vv*Tn,' theory',uu*vv/(uu**2-vv**2))
    assert abs(Tn-uu/(uu**2-vv**2))<1e-4
# 一般初期配置：A at (-a, -b) (Bの後方) s0=a
uu,vv=2.0,1.0; a,b=0.6,0.8
def f(t,y):
    xa,ya=y; dx,dy=vv*t-xa,-ya; r=np.hypot(dx,dy); return [uu*dx/r,uu*dy/r]
ev=lambda t,y: np.hypot(vv*t-y[0],y[1])-1e-6; ev.terminal=True
s=solve_ivp(f,[0,50],[-a,-b],events=ev,rtol=1e-11,atol=1e-13,max_step=1e-3)
print('general T',s.t_events[0][0],' theory',(uu*1.0+vv*a)/(uu**2-vv**2))
assert abs(s.t_events[0][0]-(uu+vv*a)/(uu**2-vv**2))<1e-4
# u=v の極限
uu=vv=1.0
s=run(uu,vv,tmax=200)
xa,ya=s.y[:,-1]; tt=s.t[-1]
print('u=v: final distance',np.hypot(vv*tt-xa,ya),' theory 0.5')
assert abs(np.hypot(vv*tt-xa,ya)-0.5)<1e-3
# 追いつく直前の向き（A→Bがx軸となす角）
s=run(2.0,1.0)
xa,ya=s.y[:,-2]; tt=s.t[-2]
print('angle near catch (deg)',np.degrees(np.arctan2(-ya,1.0*tt-xa)))
