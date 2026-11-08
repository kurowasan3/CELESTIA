# J019-1 検証：一定の力で引き上げる箱の中で糸を切った小球（相対加速度 F/M と重心の運動）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
M,m,F,g,h,u,t=sp.symbols('M m F g h u t',positive=True)
aBox=F/M-g; aBall=-g
arel=sp.simplify(aBox-aBall); print('relative accel (box up wrt ball):',arel)
t1=sp.sqrt(2*h/arel); print('t1 =',t1)
# 切った時刻の共通の速さ u（上向き）。付着直後の速さ = 重心速度
vG=u+(F/(M+m)-g)*t1
vbox=u+aBox*t1; vball=u-g*t1
print('momentum check:',sp.simplify((M*vbox+m*vball)/(M+m)-vG))
loss=sp.simplify(sp.Rational(1,2)*M*vbox**2+sp.Rational(1,2)*m*vball**2-sp.Rational(1,2)*(M+m)*vG**2)
print('KE loss =',sp.factor(loss))
# 切る前の加速度
print('accel before cut',F/(M+m)-g)
# 付着後に箱の底が受ける力（衝突後，一体で加速度 F/(M+m)-g）：小球が底から受ける力
print('normal on ball after sticking:',sp.simplify(m*(F/(M+m)-g)+m*g))
# 数値シミュレーション
Mn,mn,Fn,gn,hn,un=3.0,1.0,60.0,9.8,0.5,2.0
def rhs(t,y): return [y[2],y[3],Fn/Mn-gn,-gn]  # yb, yball, vb, vball
ev=lambda t,y: y[1]-(y[0]+0.0); ev.terminal=True; ev.direction=-1
s=solve_ivp(rhs,[0,5],[0,hn,un,un],events=ev,rtol=1e-12,atol=1e-12)
te=s.t_events[0][0]; yb,yl,vb,vl=s.y_events[0][0]
print('numeric t1',te,' formula',np.sqrt(2*hn*Mn/Fn))
vg=(Mn*vb+mn*vl)/(Mn+mn); print('numeric v after stick',vg,' formula',un+(Fn/(Mn+mn)-gn)*te)
print('numeric loss',0.5*Mn*vb**2+0.5*mn*vl**2-0.5*(Mn+mn)*vg**2,' formula m h F/(M+m)=',mn*hn*Fn/(Mn+mn))
