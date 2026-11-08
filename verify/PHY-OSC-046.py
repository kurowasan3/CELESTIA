# Q078-1: 自然長が無視できるばね（張力=k×長さ）につながれたおもりの鉛直面内の運動
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
k,m,g,a,t,v0=sp.symbols('k m g a t v_0',positive=True)
w=sp.sqrt(k/m); c=g/w**2
# 静止から放す：x = a cos wt, y = -c(1-cos wt)
x=a*sp.cos(w*t); y=-c*(1-sp.cos(w*t))
assert sp.simplify(m*sp.diff(x,t,2)+k*x)==0 and sp.simplify(m*sp.diff(y,t,2)+k*y+m*g)==0
# 直線：y = -c(1 - x/a)
assert sp.simplify(y+c*(1-x/a))==0
tb=sp.pi/(2*w); print('below O at t=',tb,' y=',sp.simplify(y.subs(t,tb)))
spd=sp.sqrt(sp.diff(x,t)**2+sp.diff(y,t)**2).subs(t,tb); print('speed=',sp.simplify(spd))
# 円運動：C を中心
R=sp.sqrt(a**2+c**2)
print('circle speed =',sp.simplify(w*R),' max length =',sp.simplify(c+R))
# 鉛直下向きに v0 で投げ下ろす
x2=a*sp.cos(w*t); y2=-c*(1-sp.cos(w*t))-v0/w*sp.sin(w*t)
assert sp.simplify(m*sp.diff(y2,t,2)+k*y2+m*g)==0 and sp.diff(y2,t).subs(t,0)==-v0
print('thrown down: below O at t=pi/(2w), y=',sp.simplify(y2.subs(t,tb)))
# 数値確認（2次元の運動方程式を直接積分）
kv,mv,gv,av=20.0,0.5,9.8,0.3
def f(tt,s): return [s[2],s[3],-kv/mv*s[0],-kv/mv*s[1]-gv]
wv=np.sqrt(kv/mv); cv=gv/wv**2
# (a) 静止から
s=solve_ivp(f,[0,np.pi/(2*wv)],[av,0,0,0],rtol=1e-12,atol=1e-14)
print('num at tb:',s.y[:,-1],' th y=',-cv)
assert abs(s.y[0,-1])<1e-8 and abs(s.y[1,-1]+cv)<1e-8
# (b) 円：速度を CP に垂直に（反時計回り）
Rv=np.hypot(av,cv); P=np.array([av,0.0]); C=np.array([0,-cv]); r=P-C
vel=wv*np.array([-r[1],r[0]])
s=solve_ivp(f,[0,2*np.pi/wv],[av,0,vel[0],vel[1]],rtol=1e-12,atol=1e-14,dense_output=True)
ts=np.linspace(0,2*np.pi/wv,200); d=np.hypot(s.sol(ts)[0]-C[0],s.sol(ts)[1]-C[1])
print('circle radius spread',d.min(),d.max(),Rv); assert d.max()-d.min()<1e-8
L=np.hypot(s.sol(ts)[0],s.sol(ts)[1]); print('max length',L.max(),cv+Rv)
# (c) 投げ下ろし
s=solve_ivp(f,[0,np.pi/(2*wv)],[av,0,0,-1.5],rtol=1e-12,atol=1e-14)
print('thrown: ',s.y[:,-1],' th y=',-(cv+1.5/wv)); assert abs(s.y[1,-1]+cv+1.5/wv)<1e-8
print('OK')
