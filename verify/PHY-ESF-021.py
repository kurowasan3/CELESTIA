# Q183-1 双極子のまわりの円形の壁に沿う荷電小球：垂直抗力が一定
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
th,th0,k,Q,p,m,R=sp.symbols('theta theta0 k Q p m R',positive=True)
V=lambda r,t: k*p*sp.cos(t)/r**2
r=sp.symbols('r',positive=True)
Er=-sp.diff(V(r,th),r); Eth=-sp.diff(V(r,th),th)/r
print('Er',sp.simplify(Er),' Eth',sp.simplify(Eth))
v2=2*Q*(V(R,th0)-V(R,th))/m
N=m*v2/R+Q*Er.subs(r,R)        # 内向き正: m v^2/R = N - Q E_r
print('N=',sp.simplify(N))
assert sp.simplify(N-2*k*Q*p*sp.cos(th0)/R**3)==0
print('vmax(th0=pi/2)=',sp.sqrt(sp.simplify(v2.subs({th0:sp.pi/2,th:sp.pi}))))
# 数値：有限の d をもつ2点電荷の正確な場で、壁なしの自由粒子を th0=pi/2 から放す
kk,QQ,mm,RR,d=1.0,1.0,1.0,1.0,1e-3
q0=1.0/d  # p=1
def acc(x,y):
    ax=ay=0
    for qq,xx in ((q0,d/2),(-q0,-d/2)):
        dx,dy=x-xx,y; r3=(dx*dx+dy*dy)**1.5
        ax+=kk*QQ*qq*dx/r3; ay+=kk*QQ*qq*dy/r3
    return ax/mm,ay/mm
def f(t,s):
    ax,ay=acc(s[0],s[1]); return [s[2],s[3],ax,ay]
ev=lambda t,s: s[3] if t>0.1 else 1; 
s=solve_ivp(f,[0,6],[0,1,0,0],rtol=1e-11,atol=1e-12,dense_output=True)
rr=np.hypot(s.y[0],s.y[1]); print('free particle r range',rr.min(),rr.max())
thn=np.arctan2(s.y[1],s.y[0])%(2*np.pi); print('theta range',thn.min(),thn.max())
sp_=np.hypot(s.y[2],s.y[3]); print('vmax num',sp_.max(),' theory',np.sqrt(2))
# th0=pi/3 では壁の力が一定 kQp/R^3 になることを、壁拘束の数値（角度だけの運動方程式）で確認
def g(t,s):
    t_,w=s; return [w, (QQ*kk*1*np.sin(t_)/RR**3)/mm/RR]
s2=solve_ivp(g,[0,5],[np.pi/3,0],rtol=1e-11,atol=1e-12)
Nnum=mm*RR*s2.y[1]**2+QQ*2*kk*np.cos(s2.y[0])/RR**3
print('N range (th0=60deg)',Nnum.min(),Nnum.max(),' theory',kk*QQ/RR**3)
print('theta max',s2.y[0].max(),' theory',2*np.pi-np.pi/3)
print('OK')
