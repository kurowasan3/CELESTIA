# S0115-1 赤道で真上に投げ上げた物体の着地点のずれ（面積速度一定）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
t,v0,g,W=sp.symbols('t v0 g Omega',positive=True)
y=v0*t-g*t**2/2; Tf=v0/g
# 地面に沿った東向き相対変位 = ∫ -2 W y dt
up=sp.integrate(-2*W*y,(t,0,Tf)); tot=sp.integrate(-2*W*y,(t,0,2*Tf))
print('top offset',sp.simplify(up),' landing',sp.simplify(tot))
h=v0**2/(2*g); print('rel vel at top',sp.simplify(-2*W*h))
# 頂点から静かに落とした場合（比較）: ∫ 2W(h-z) dt, z=h-g t^2/2
drop=sp.integrate(W*g*t**2,(t,0,Tf)); print('drop east',sp.simplify(drop))
# 独立検証: 慣性系で厳密な万有引力のもと数値積分
R=6.4e6; Om=7.3e-5; gv=9.8; g0=gv+Om**2*R   # 地表の見かけの重力 = g0 - Om^2 R
GM=g0*R**2
def f(tt,Y):
    x,yy,vx,vy=Y; rr=np.hypot(x,yy); a=-GM/rr**3
    return [vx,vy,a*x,a*yy]
def land(tt,Y): return np.hypot(Y[0],Y[1])-R
land.terminal=True; land.direction=-1
for vv in [49.0,98.0]:
    # 地点 (R,0)，鉛直上向き vv，東向き(=+y方向，反時計回り) Om R
    s=solve_ivp(f,[1e-9,4*vv/gv],[R,0,vv,Om*R],events=land,rtol=1e-13,atol=1e-9,first_step=1e-3)
    tl=s.t_events[0][0]; xl,yl=s.y_events[0][0][:2]
    ang=np.arctan2(yl,xl); d=R*(ang-Om*tl)
    th=-(4/3)*Om*vv**3/gv**2
    print('v0',vv,'numeric east offset',d,' theory',th)
    assert abs(d-th)/abs(th)<0.02

# ---- (6) latitude lambda: west offset (4/3) W v0^3 cos(lam)/g^2 ; 3D inertial integration ----
def f3(tt,Y):
    x,yy,z,vx,vy,vz=Y; rr=np.sqrt(x*x+yy*yy+z*z); a=-GM/rr**3
    return [vx,vy,vz,a*x,a*yy,a*z]
def land3(tt,Y): return np.sqrt(Y[0]**2+Y[1]**2+Y[2]**2)-R
land3.terminal=True; land3.direction=-1
for latd in [0.0,40.0,75.0,90.0]:
    lat=np.radians(latd); vv=98.0
    P=np.array([R*np.cos(lat),0,R*np.sin(lat)]); up=P/R
    east=np.array([0,1.0,0])
    # apparent vertical differs slightly from radial by centrifugal; use radial 'up' (first-order effect on E-W negligible)
    V0=vv*up+Om*R*np.cos(lat)*east
    s=solve_ivp(f3,[1e-9,4*vv/gv],list(P)+list(V0),events=land3,rtol=1e-13,atol=1e-9,first_step=1e-3)
    tl=s.t_events[0][0]; X=s.y_events[0][0][:3]
    lon=np.arctan2(X[1],X[0])-Om*tl      # longitude relative to rotating P
    dE=R*np.cos(lat)*lon
    th=-(4/3)*Om*vv**3*np.cos(lat)/gv**2
    print('(6) lat %4.0f numeric east offset %.4f m theory %.4f m'%(latd,dE,th))
    assert abs(dE-th)<0.02*abs(-(4/3)*Om*vv**3/gv**2)
