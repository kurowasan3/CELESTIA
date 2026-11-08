# J078-1 検証：点波源からの屈折波の波面（境界上の交点の速さと見かけの波源）
import sympy as sp, numpy as np
from scipy.optimize import minimize_scalar
V,d,t,a=sp.symbols('V d t alpha',positive=True)
# (1) x_T(t), u_T
xT=sp.sqrt(V**2*t**2-d**2); uT=sp.diff(xT,t)
sini=xT/(V*t); print('uT*sin i =',sp.simplify(uT*sini))
# (2) sin theta = (V/2)/uT
print('sin theta / sin i =',sp.simplify((V/2/uT)/sini))
# 波面のパラメータ表示（入射角 alpha の射線）
b=sp.asin(sp.sin(a)/2)
tau=t-d/(V*sp.cos(a)); s=V/2*tau
X=d*sp.tan(a)+s*sp.sin(b); Y=-s*sp.cos(b)
t1=2*d/V
Xs=X.subs(t,t1); Ys=Y.subs(t,t1)
print('S:',sp.simplify(Xs.subs(a,0)),sp.simplify(Ys.subs(a,0)))
# T: s=0 -> cos a = d/(V t1) = 1/2 -> a=pi/3
aT=sp.pi/3; print('T:',sp.simplify(Xs.subs(a,aT)),sp.simplify(Ys.subs(a,aT)),' sin theta:',sp.simplify(sp.sin(b).subs(a,aT)))
# (4) 曲率半径 at alpha=0
dX=sp.diff(Xs,a); dY=sp.diff(Ys,a); ddX=sp.diff(Xs,a,2); ddY=sp.diff(Ys,a,2)
kap=sp.simplify(((dX*ddY-dY*ddX)/(dX**2+dY**2)**sp.Rational(3,2)).subs(a,0))
print('curvature radius at S:',sp.simplify(1/kap))
# 射線の逆延長の y 切片（一般の alpha）
yint=d*sp.tan(a)/sp.tan(b)
print('apparent y (alpha->0):',sp.limit(yint,a,0),'  at T:',sp.nsimplify(sp.simplify(yint.subs(a,aT))))
print('QT =',sp.sqrt(3*d**2+4*d**2),' vs R=5d/2')
# 数値：フェルマー（最短時間）で波面を求めパラメータ表示と比較
Vn,dn=1.0,1.0; tn=2.0
def tarr(x,y):
    f=lambda xb: np.hypot(xb,dn)/Vn+np.hypot(x-xb,y)/(Vn/2)
    r=minimize_scalar(f,bounds=(0,max(x,1e-9)) if x>0 else (-1,1),method='bounded',options={'xatol':1e-12}); return r.fun
Xf=sp.lambdify(a,Xs.subs({V:1,d:1})); Yf=sp.lambdify(a,Ys.subs({V:1,d:1}))
for al in [0.0,0.3,0.6,0.9,1.0]:
    x,y=Xf(al),Yf(al); print(f'alpha={al}: point ({x:.4f},{y:.4f}) arrival time {tarr(x,y):.6f} (should be 2)')
