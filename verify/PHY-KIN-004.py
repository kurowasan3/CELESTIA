# T20141-1 : 水平の長さが一定のジャンプ台（θ で発射の速さが変わる）と最大到達距離
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
g,H,b,th,t=sp.symbols('g H b theta t',positive=True)
v2=2*g*(H-b*sp.tan(th))
s=sp.simplify(v2*sp.sin(2*th)/g)
print('s(theta)=',s)
# t=tanθ に変数変換して最大
T=sp.symbols('T',positive=True)
sT=4*(H-b*T)*T/(1+T**2)
crit=sp.solve(sp.diff(sT,T),T); print('tan theta* =',crit)
Ts=[c for c in crit if c.is_positive!=False][0]
print('tan2θ* =',sp.simplify(2*Ts/(1-Ts**2)))
smax=sp.simplify(sT.subs(T,Ts)); print('s_max=',smax,' = 2H tanθ* ?',sp.simplify(smax-2*H*Ts))
print('H=sqrt3 b: theta*=',sp.simplify(sp.atan(Ts.subs(H,sp.sqrt(3)*b))),' smax=',sp.simplify(smax.subs(H,sp.sqrt(3)*b)))
# H>>b での 45°-θ*
eps=sp.symbols('epsilon',positive=True)
thstar=sp.atan(Ts).subs(b,eps*H)
print('series 45°-θ* :',sp.series(sp.pi/4-thstar,eps,0,3))
# 数値: 実際に放物運動を積分して最大到達距離と角度を確認 (H=√3, b=1, g=9.8)
Hv,bv,gv=np.sqrt(3),1.0,9.8
def rng(theta):
    vv=np.sqrt(2*gv*(Hv-bv*np.tan(theta)))
    ev=lambda t,s:s[1]; ev.terminal=True; ev.direction=-1
    sol=solve_ivp(lambda t,s:[s[2],s[3],0,-gv],[0,10],[0,0,vv*np.cos(theta),vv*np.sin(theta)],events=ev,rtol=1e-11,atol=1e-12,first_step=1e-4)
    return sol.y_events[0][0][0]
r=minimize_scalar(lambda th:-rng(th),bounds=(0.1,1.0),method='bounded')
print('numeric optimum deg',np.degrees(r.x),' smax',-r.fun,' expected 30, 2b=2')
# 条件 H>b tanθ* が自動的に成り立つ: H-b tanθ* = H(1+tan²)/2
print('H-b tan* =',sp.simplify(H-b*Ts-H*(1+Ts**2)/2))
