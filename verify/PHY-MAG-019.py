# O20072-1 : 鉛直磁場中の斜面上の荷電小物体（ドリフトと面から離れる点）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
q,m,B,g,th,u,t=sp.symbols('q m B g theta u t',positive=True)
vx,vy=sp.symbols('v_x v_y')
# 斜面座標: x 水平, y 斜面上向き, z 法線. 鉛直上向き = (0, sin, cos)
up=sp.Matrix([0,sp.sin(th),sp.cos(th)]); Bv=-B*up
v=sp.Matrix([vx,vy,0]); FL=q*v.cross(Bv)
print('Lorentz=',sp.simplify(FL.T))
N=sp.simplify(m*g*sp.cos(th)-FL[2]); print('N=',N)
ax=FL[0]/m; ay=(FL[1]-m*g*sp.sin(th))/m
vd=sp.solve(ay.subs(vy,0),vx)[0]; print('drift vd=',sp.simplify(vd))
w=sp.symbols('w'); om=q*B*sp.cos(th)/m
VX=vd+w*sp.cos(om*t); VY=w*sp.sin(om*t)
print('check eqs:',sp.simplify(sp.diff(VX,t)-ax.subs({vx:VX,vy:VY})),sp.simplify(sp.diff(VY,t)-ay.subs({vx:VX,vy:VY})))
cc=sp.solve(N.subs(vx,-sp.Symbol('C')),sp.Symbol('C'))[0]; print('c=',sp.simplify(cc))
umax=sp.simplify(2*vd+cc); print('u_max=',umax, '=', sp.simplify(umax - m*g/(q*B)*(2*sp.tan(th)+1/sp.tan(th))))
# 離れる点
wv=u-vd; cosT=(-cc-vd)/wv
Y=sp.simplify(wv/om*(1-cosT)); print('Y_leave=',Y,' check (u+c)/om:',sp.simplify(Y-(u+cc)/om))
v2=sp.simplify(cc**2+wv**2-(cc+vd)**2); E=sp.simplify(u**2-2*g*sp.sin(th)*Y)
print('speed^2 kin =',sp.factor(v2),' energy =',sp.factor(E),' diff',sp.simplify(v2-E))
# 数値: 運動方程式を積分し N=0 の点を探す
qv,mv,Bv_,gv,thv=1.0,1.0,1.0,9.8,np.pi/6
vdn=float(vd.subs({q:qv,m:mv,B:Bv_,g:gv,th:thv})); cn=float(cc.subs({q:qv,m:mv,B:Bv_,g:gv,th:thv}))
for uv in [2*vdn+cn-1, 2*vdn+cn+3]:
    def f(tt,y):
        X,Yy,VX_,VY_=y
        return [VX_,VY_,-qv*Bv_*np.cos(thv)/mv*VY_, qv*Bv_*np.cos(thv)/mv*VX_-gv*np.sin(thv)]
    def ev(tt,y): return mv*gv*np.cos(thv)+qv*Bv_*np.sin(thv)*y[2]
    ev.terminal=True; ev.direction=-1
    s=solve_ivp(f,(0,20),[0,0,uv,0],events=ev,rtol=1e-11,atol=1e-12)
    if s.t_events[0].size:
        yl=s.y_events[0][0]; om_=qv*Bv_*np.cos(thv)/mv
        print('u=',uv,' leaves at y=',yl[1],' theory',(uv+cn)/om_,' speed',np.hypot(yl[2],yl[3]),' theory',np.sqrt(uv**2-2*gv*np.sin(thv)*(uv+cn)/om_))
    else: print('u=',uv,' never leaves (theory: u<=umax)', uv<=2*vdn+cn, ' max y',max(s.y[1]),' theory 2w/om',2*(uv-vdn)/(qv*Bv_*np.cos(thv)/mv))
