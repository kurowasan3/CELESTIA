# P14-1: 上部を切り取ったループの飛び越し（切り口の角度と最低出発高さ）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
a,R,g,h,s=sp.symbols('alpha R g h s',positive=True)
# P=(R cos a, R sin a) (O原点), 速度向き (-sin a, cos a)
v2=sp.symbols('v2',positive=True)
t=sp.symbols('t',positive=True)
vx=-sp.sqrt(v2)*sp.sin(a); vy=sp.sqrt(v2)*sp.cos(a)
tf=2*vy/g
sol=sp.solve(sp.Eq(-vx*tf,2*R*sp.cos(a)),v2); print('v_P^2 =',sp.simplify(sol[0]))
hP=sp.simplify(R+R*sp.sin(a)+sol[0]/(2*g)); print('h =',hP)
da=sp.solve(sp.diff(hP,a),a); print('opt alpha',da, 'hmin',sp.simplify(hP.subs(a,sp.pi/4)))
# N>=0 条件: 2(h-R)/R >= 3 sin a
print('N at P margin/(mg):',sp.simplify((2*(hP-R)/R-3*sp.sin(a))))
# 2根の積
H=sp.symbols('H',positive=True)
print(sp.Poly(sp.expand((s+1/(2*s)-(H-1))*2*s),s).all_coeffs())
# 数値: alpha=45度, h=(1+sqrt2)R で P から放物運動し Q に到達, Q での速度が接線方向か
Rn,gn=1.0,9.8; al=np.pi/4; hn=(1+np.sqrt(2))*Rn
vP=np.sqrt(2*gn*(hn-Rn-Rn*np.sin(al)))
y0=[Rn*np.cos(al),Rn*np.sin(al),-vP*np.sin(al),vP*np.cos(al)]
# 放物線は Q で円に接する（符号変化しない）ので，x=-R cos a を通過する時刻で判定
ev=lambda t,y: y[0]+Rn*np.cos(al)
ev.terminal=True
r=solve_ivp(lambda t,y:[y[2],y[3],0,-gn],[0,5],y0,events=ev,rtol=1e-11,atol=1e-12)
x,y,u,w=r.y[:,-1]; print('landing',x,y,'Q',-Rn*np.cos(al),Rn*np.sin(al),'radial vel',(x*u+y*w)/Rn)
ts=np.linspace(0,r.t[-1],2001); yy=r.sol(ts) if r.sol is not None else None
# 飛行中は円の外側（切り取った部分）を通り，残った軌道と交わらないことの確認
for al in [np.pi/6,np.pi/4,np.pi/3,1.3]:
    hn=Rn*(1+np.sin(al)+1/(2*np.sin(al))); vP=np.sqrt(2*gn*(hn-Rn-Rn*np.sin(al)))
    T=2*vP*np.cos(al)/gn; ts=np.linspace(1e-6,T-1e-6,4001)
    X=Rn*np.cos(al)-vP*np.sin(al)*ts; Y=Rn*np.sin(al)+vP*np.cos(al)*ts-gn*ts**2/2
    rr=np.hypot(X,Y); ang=np.arctan2(Y,np.abs(X))
    print('alpha',al,'min r/R',rr.min(),'min elevation',ang.min(),'>=alpha',ang.min()>=al-1e-9)
# (5) 2根がともに (0,1) にある範囲の確認
for Hh in [2.42,2.45,2.49,2.5,2.6]:
    rts=np.roots([2,-2*(Hh-1),1]); print('h/R',Hh,'roots',np.round(rts,4),'prod',np.prod(rts).real)
