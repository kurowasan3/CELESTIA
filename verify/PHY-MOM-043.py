# PHY-MOM-043 検証: モーターで回るおもりと振り子の2回衝突
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
e,u,g,l,mp=sp.symbols('e u g ell mprime',positive=True)
# 衝突: おもりの速度は u のまま。反発係数はおもりの系で適用
v1=sp.symbols('v1'); v1s=sp.solve(sp.Eq(v1-u,-e*(0-u)),v1)[0]
v2=sp.symbols('v2'); v2s=sp.solve(sp.Eq(v2-u,-e*(-v1s-u)),v2)[0]
print('v1=',sp.factor(v1s),' v2=',sp.factor(v2s))
us=sp.solve(sp.Eq(v1s**2,2*g*l),u)[0]; print('u=',us)
sb=sp.symbols('s')  # sin(beta)
sbs=sp.solve(sp.Eq(v2s.subs(u,us)**2-2*g*l*(1+sb),g*l*sb),sb)[0]; print('sin beta=',sp.factor(sbs))
print('e_max=',sp.solve(sp.Eq(sbs,1),e),' e(beta=30deg)=',sp.solve(sp.Eq(sbs,sp.Rational(1,2)),e))
J=mp*(v2s+v1s); W=sp.simplify(J*u); dK=sp.Rational(1,2)*mp*(v2s**2-v1s**2)
print('W=',sp.factor(W.subs(u,us)),' loss=',sp.factor(sp.simplify((W-dK).subs(u,us))))
# 数値: 振り子の運動を積分し，2回目の後に張力0となる角度を確認
for ev in [0.2,0.323,0.5]:
    gv,lv=9.8,1.0; uv=float(us.subs({g:gv,l:lv,e:ev})); v=(1+ev)**2*uv
    def rhs(t,y): return [y[1], -gv/lv*np.sin(y[0])]
    def tension(t,y): return lv*y[1]**2+gv*np.cos(y[0])   # T/m
    tension.terminal=True; tension.direction=-1
    s=solve_ivp(rhs,[0,5],[0,v/lv],events=tension,rtol=1e-11,atol=1e-12)
    th=s.y_events[0][0][0]; beta=th-np.pi/2
    print(ev,'num sin(beta)=',np.sin(beta),' formula',2*ev*(2+ev)/3)
    # 1回目: 水平まで上がり張力が正のまま
    s1=solve_ivp(rhs,[0,5],[0,(1+ev)*uv/lv],events=tension,rtol=1e-11,atol=1e-12,max_step=1e-3)
    print('   first swing max angle(deg)',np.degrees(max(s1.y[0])),' slack events',len(s1.t_events[0]))
