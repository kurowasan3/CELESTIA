# J001-1 検証：流速 v が船の速さ u より大きい川の横断（最小の流され距離と2つの船首方向）
import sympy as sp, numpy as np
u,v,W,b=sp.symbols('u v W beta',positive=True)
# beta: 船首と「上流向きの岸」のなす角．x 下流向き，y 対岸向き
vx=v-u*sp.cos(b); vy=u*sp.sin(b)
X=W*vx/vy
# (1) 最小の流され距離
crit=sp.solve(sp.diff(X,b),b); print('crit beta:',crit)
bc=sp.acos(u/v)
Xmin=sp.simplify(X.subs(b,bc)); print('Xmin =',Xmin)
tmin=sp.simplify((W/vy).subs(b,bc)); print('time at Xmin =',tmin)
# 合成速度と船首の直交
print('dot(vel, heading) at bc =',sp.simplify((vx*(-sp.cos(b))+vy*sp.sin(b)).subs(b,bc)))
# (3) v=2u, L=2W
eq=sp.Eq((X/W).subs(v,2*u),2)
t=sp.symbols('t',positive=True)
sols=sp.solve(sp.Eq((2-(1-t**2)/(1+t**2))/((2*t)/(1+t**2)),2),t); print('tan(b/2)=',sols)
for s in sols:
    bb=2*sp.atan(s); print(' beta=',sp.N(bb*180/sp.pi),' sin=',sp.simplify(sp.sin(bb)),' cos=',sp.simplify(sp.cos(bb)),' time=',sp.simplify(W/(u*sp.sin(bb))))
# (4) 1方向しかない => L=Xmin
uu=sp.solve(sp.Eq(Xmin,2*W),u); print('u for unique heading:',uu)
# 数値で：beta を細かく走査し最小値
un,vn,Wn=1.0,2.0,1.0
bs=np.linspace(0.01,np.pi-0.01,200001); Xs=Wn*(vn-un*np.cos(bs))/(un*np.sin(bs))
print('numeric Xmin',Xs.min(),' formula',np.sqrt(3))
# 上流側に着けないこと：vx>0 for u<v
print('vx min over beta =',(vn-un*np.cos(bs)).min())
