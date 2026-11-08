# Q074-1: 樋に当たる水流の分かれ方と台車の運動
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
mu,v,al,M0,t,u=sp.symbols('mu v alpha M_0 t u',positive=True)
m1,m2=sp.symbols('m1 m2')
sol=sp.solve([m1+m2-mu, m1*v-m2*v-mu*v*sp.sin(al)],[m1,m2])
print(sol); assert sp.simplify(sol[m1]/sol[m2]-(1+sp.sin(al))/(1-sp.sin(al)))==0
# 手で支える力：バケツに入る水平運動量の流量
Fh=sol[m1]*v-sol[m2]*v; assert sp.simplify(Fh-mu*v*sp.sin(al))==0
# 台車：(M0+mu t)u = mu v sin(al) t
U=mu*v*sp.sin(al)*t/(M0+mu*t)
lhs=sp.diff((M0+mu*t)*U,t); assert sp.simplify(lhs-mu*v*sp.sin(al))==0
# 準定常の別解：台車系で水平運動量の流入 mu (v sin(al) - u) が (M0+mu t) du/dt に等しい（質量増加分は相対速度0で取り込まれる…ではなく相対速度で入る）
# d/dt[(M0+mu t)u] = mu v sin(al)  <=>  (M0+mu t) du/dt = mu (v sin(al) - u)
ode=lambda tt,y,p: [p['mu']*(p['v']*np.sin(p['al'])-y[0])/(p['M0']+p['mu']*tt)]
p=dict(mu=2.0,v=5.0,al=np.pi/6,M0=40.0)
s=solve_ivp(ode,[0,100],[0.0],args=(p,),rtol=1e-10,atol=1e-12,dense_output=True)
for tt in [5,20,100]:
    th=float(U.subs({mu:2,v:5,al:sp.pi/6,M0:40,t:tt})); print(tt,s.sol(tt)[0],th); assert abs(s.sol(tt)[0]-th)<1e-6
# 台車の速さ u のときの流量の差（台車系）
w=sp.sqrt((v*sp.sin(al)-u)**2+(v*sp.cos(al))**2)
diff=mu*(v*sp.sin(al)-u)/w
val=sp.simplify(diff.subs({al:sp.pi/6,u:v/4}))
print('u=v/4, alpha=30deg: (m1-m2)/mu =',sp.nsimplify(val/mu), float(val/mu))
r=sp.simplify(((mu+diff)/(mu-diff)).subs({al:sp.pi/6,u:v/4}))
print('ratio =',sp.nsimplify(sp.radsimp(r)), float(r))
# 時刻：u=v sin(al)/2 となる時刻 t = M0/mu
print(sp.solve(sp.Eq(U,v*sp.sin(al)/2),t))
print('OK')
