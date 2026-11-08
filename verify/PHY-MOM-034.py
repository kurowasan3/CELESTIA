# Q103-1 穴から落ちる鎖（山からの取り込み・床での停止・自由落下の3段階）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,h,x,s,t=sp.symbols('g h x s t',positive=True)
# 段階1: d(xv)/dt = x g -> v^2=2gx/3 を確認
v2=2*g*x/3; v=sp.sqrt(v2)
# d(xv)/dt = v d(xv)/dx
print('phase1 check',sp.simplify(v*sp.diff(x*v,x)-x*g))
t1=sp.sqrt(6*h/g); print('t1 check x(t1)=',sp.simplify((g/6)*t1**2))
# エネルギー損失の割合
PEdrop=g*x**2/2; KE=x*v2/2
print('lost fraction',sp.simplify((PEdrop-KE)/PEdrop))
# 段階2: d(v^2)/ds = 2g - 2 v^2/h, v^2(0)=2gh/3
Y=sp.Function('Y')
sol=sp.dsolve(sp.Eq(Y(s).diff(s),2*g-2*Y(s)/h),Y(s),ics={Y(0):2*g*h/3})
print('phase2',sp.simplify(sol.rhs))
vend2=sp.simplify(sol.rhs.subs(s,h)); print('v^2 at rear leaving',vend2)
Ffloor=sp.simplify(vend2+g*h); print('floor force/lambda',sp.factor(Ffloor))
# 段階3
vv=sp.sqrt(vend2); t3=(-vv+sp.sqrt(vv**2+2*g*h))/g
print('t3',sp.simplify(t3), float(t3.subs({g:9.8,h:1})))
# ---- 数値: 離散リンクモデル（質量 dm の玉が順に取り込まれる） ----
gv,hv=9.8,1.0; N=4000; L=2*hv; dl=L/N
# 状態: 吊り下がり部分の長さ X(<=h), 速度 V, 通過した長さ S
# 取り込み: 運動量保存（完全非弾性） / 床到達: 下端の玉は床で止まる（吊り部分の他の玉に影響なし）
X=dl; V=0.0; S=dl; T=0.0
dt=2e-5; floor_hit=None; rear=None
while S< L-1e-12:
    # 吊り部分は重力で自由に加速（内部張力は下端の玉を止めるのに使われない）
    V+=gv*dt; T+=dt
    move=V*dt
    # 取り込み: move だけ新たに通過 -> 運動量保存
    newlen=min(move,L-S)
    mass_h=min(X,hv)
    V=V*mass_h/(mass_h+newlen)
    X+=newlen; S+=newlen
    if X>hv:
        if floor_hit is None: floor_hit=(T,V)
        X=hv
print('front hits floor t,v sim',floor_hit,' theory',float(t1.subs({g:gv,h:hv})),np.sqrt(2*gv*hv/3))
print('rear leaves v sim',V,' theory',float(sp.sqrt(vend2).subs({g:gv,h:hv})))
# 段階1を ODE で: x0 を小さくして t1 が sqrt(6h/g) に収束することを確認
for x0 in [1e-3,1e-5,1e-7]:
    ev=lambda t,y: y[0]-hv; ev.terminal=True
    so=solve_ivp(lambda t,y:[y[1], gv - y[1]**2/y[0]],[0,5],[x0,0.0],events=ev,rtol=1e-11,atol=1e-14)
    print('x0',x0,'t1',so.t_events[0][0],'v',so.y_events[0][0][1])
# 段階2を ODE で
so2=solve_ivp(lambda t,y:[y[1], gv - y[1]**2/hv],[0,5],[0,np.sqrt(2*gv*hv/3)],events=(lambda t,y:y[0]-hv),rtol=1e-11,atol=1e-13)
print('phase2 v at s=h',so2.y_events[0][0][1])
