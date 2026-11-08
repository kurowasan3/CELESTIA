# P13 類題：静止を保つ回転枠の板の上を走る小動物（力の作用線が軸を通る→単振動／不安定）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,g,d,a,x,t,M,v0,b=sp.symbols('m g d a x t M v0 b',positive=True)
lam=sp.symbols('lam')
# 板が軸の下 d：動物の位置 (x,-d)、Oは原点。力 F = lam*(O - P) = lam*(-x, d)
Fx=-lam*x; Fy=lam*d
lsol=sp.solve(sp.Eq(Fy,m*g),lam)[0]
print('ax =',sp.simplify(Fx.subs(lam,lsol)/m))
# トルク（Oまわり）: 動物が板に及ぼす力 -F を点 (x,-d) に
tq=sp.simplify(x*(-Fy) - (-d)*(-Fx)).subs(lam,lsol); print('torque on frame',sp.simplify(tq))
w=sp.sqrt(g/d)
xt=a*sp.cos(w*t); assert sp.simplify(sp.diff(xt,t,2)+g/d*xt)==0
print('friction ratio |Fx|/Fy max = a/d')
# 板が軸の上 d：F = lam*(-x,-d), Fy=mg -> lam=-mg/d, ax = +g x/d
lam2=sp.solve(sp.Eq(-lam*d,m*g),lam)[0]; print('above: ax=',sp.simplify(-lam2*x/m))
W=sp.sqrt(g/d)
xs=v0/(2*W)*(sp.exp(W*t)-sp.exp(-W*t))
assert sp.simplify(sp.diff(xs,t,2)-W**2*xs)==0
T=sp.log(b*W/v0+sp.sqrt(1+(b*W/v0)**2))/W
print('time to reach b:',T)
print('check x(T)=b:',sp.simplify(xs.subs(t,T)-b))
# numerical: d=0.1, g=9.8, a=0.03: animal motion with force through O; check period
dn,gn,an=0.10,9.8,0.03
f=lambda t,y:[y[1],-gn*y[0]/dn]
s=solve_ivp(f,[0,2],[an,0],max_step=1e-4,rtol=1e-10,dense_output=True)
# find first return to x=a
ts=np.linspace(0.01,2,200001); xs_=s.sol(ts)[0]; i=np.argmax((xs_[1:]<xs_[:-1]) & (np.arange(len(ts)-1)>1000) & (xs_[:-1]>an-1e-6))
print('period theory',2*np.pi*np.sqrt(dn/gn))
# above case numeric: v0=0.2, b=0.3, d=0.1
v0n,bn=0.2,0.3
f2=lambda t,y:[y[1],gn*y[0]/dn]
ev=lambda t,y:y[0]-bn; ev.terminal=True
s2=solve_ivp(f2,[0,5],[0,v0n],events=ev,rtol=1e-11,atol=1e-13)
Tth=float(T.subs({v0:v0n,b:bn,g:gn,d:dn}))
print('time to edge numeric',s2.t_events[0][0],' theory',Tth); assert abs(s2.t_events[0][0]-Tth)<1e-6
print('numbers: T=',2*np.pi*np.sqrt(0.2/9.8),' amax(mu=0.6,d=0.2)=',0.6*0.2)
