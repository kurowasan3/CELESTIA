# T20121-1 : 台地を横切る小球の「屈折」
import sympy as sp, numpy as np
from scipy.optimize import minimize_scalar
g,h,v0,i,w=sp.symbols('g h v0 i w',positive=True)
# (1) y成分保存, x成分のエネルギー
vx_H=sp.sqrt(v0**2*sp.cos(i)**2-2*g*h)
v1=sp.sqrt(v0**2-2*g*h)
assert sp.simplify(vx_H**2+(v0*sp.sin(i))**2-v1**2)==0
sinr=v0*sp.sin(i)/v1
print('sin i/sin r =',sp.simplify(sp.sin(i)/sinr))
# (4) v0=2sqrt(gh)
sub={v0:2*sp.sqrt(g*h)}
ic=sp.acos(sp.sqrt(2*g*h)/(2*sp.sqrt(g*h)))
print('critical i =',sp.simplify(ic))
s=sp.simplify(sinr.subs(sub).subs(i,sp.pi/6)); print('sin r at 30deg =',s)
r=sp.asin(s); dy=sp.simplify(w*sp.tan(r)); print('dy on H =',dy)
shift=sp.simplify(dy-w*sp.tan(sp.pi/6)); print('shift =',shift, float(shift.subs(w,1)))
# (5) 下りるときの実効屈折率 n(u)=sqrt(1+2gh/u^2)
u=sp.symbols('u',positive=True)
n=sp.sqrt(u**2+2*g*h)/u; print('n(u)=',sp.simplify(n),' dn/du<0:',sp.simplify(sp.diff(n,u)))
# (6) 極値原理の数値確認: A=(-1,0) on L, B=(1,Y) on H, boundary x=0 (台地の西端)
G,Hh,V0=9.8,0.5,4.0; V1=np.sqrt(V0**2-2*G*Hh)
A=np.array([-1.0,0.0]); B=np.array([1.0,1.3])
L1=lambda y: np.hypot(1,y); L2=lambda y: np.hypot(1,B[1]-y)
S=lambda y: V0*L1(y)+V1*L2(y); T=lambda y: L1(y)/V0+L2(y)/V1
for name,f in [('S',S),('T',T)]:
    y=minimize_scalar(f,bounds=(-5,5),method='bounded').x
    si=y/L1(y); sr=(B[1]-y)/L2(y)
    print(name,'stationary: v0 sin i - v1 sin r =',V0*si-V1*sr,' sin i/v0 - sin r/v1 =',si/V0-sr/V1)
# 数値積分: なめらかな斜面 z=f(x)（幅 e の smoothstep）上の束縛運動を解き，上面Hでの sin r を確認
from scipy.integrate import solve_ivp
G,Hh=9.8,0.5; e=0.05; V0=2*np.sqrt(G*Hh)
def f(x):
    s=np.clip(x/e,0,1); return Hh*(3*s**2-2*s**3)
def fp(x):
    s=np.clip(x/e,0,1); return Hh*(6*s-6*s**2)/e if 0<x<e else 0.0
def fpp(x):
    s=np.clip(x/e,0,1); return Hh*(6-12*s)/e**2 if 0<x<e else 0.0
def rhs(t,st):
    x,y,vx,vy=st
    a=-(fp(x)*fpp(x)*vx**2+G*fp(x))/(1+fp(x)**2)
    return [vx,vy,a,0]
for deg,exp in [(30,1/np.sqrt(2)),(20,None),(50,'reflect')]:
    ii=np.radians(deg)
    ev=lambda t,s:s[0]-0.2; ev.terminal=True
    ev2=lambda t,s:s[0]+0.1; ev2.terminal=True
    sol=solve_ivp(rhs,[0,5],[-0.05,0,V0*np.cos(ii),V0*np.sin(ii)],max_step=2e-4,events=[ev,ev2],rtol=1e-10,atol=1e-12)
    vx,vy=sol.y[2,-1],sol.y[3,-1]
    th_sinr=V0*np.sin(ii)/np.sqrt(V0**2-2*G*Hh) if V0*np.cos(ii)>np.sqrt(2*G*Hh) else None
    print(deg,'deg: final x=%.3f vx=%.4f vy=%.4f sin=%.5f theory sin r=%s'%(sol.y[0,-1],vx,vy,vy/np.hypot(vx,vy),th_sinr))
