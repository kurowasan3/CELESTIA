# PHY-WAV-024 検証：平面波と円形波の山の交点（高い山）の運動
import sympy as sp
import numpy as np
from scipy.optimize import brentq
x,y,lam,v,t,c=sp.symbols('x y lambda v t c',positive=True)
# (1) r - y = c の曲線
r=sp.sqrt(x**2+y**2)
Y=sp.solve(sp.Eq(r-y,c),y); print('(1) y =',Y)
# (3) P(4λ,3λ) の速度: 陰関数微分
xs,ys,rs=sp.symbols('xs ys rs')
# x xdot + y ydot = r rdot, ydot=v, rdot=v
xd=sp.solve(sp.Eq(4*lam*sp.Symbol('xd')+3*lam*v,5*lam*v),sp.Symbol('xd'))[0]; print('(3) xdot =',xd,' speed',sp.sqrt(xd**2+v**2))
# 数値で確認：λ=1, v=1, 時刻 t における山の交点（n=2 の曲線上，P から出発）を直接解く
L,V=1.0,1.0
def crest_point(tt):
    # 平面波の山: y = 3 + V tt（P を通る山），円形波の山: r = 5 + V tt
    yy=3+V*tt; rr=5+V*tt
    return np.sqrt(rr**2-yy**2), yy
h=1e-6
for tt in [0,1,5,12]:
    (x1,y1),(x2,y2)=crest_point(tt-h),crest_point(tt+h)
    print(' t=',tt,' pos',crest_point(tt),' vel',((x2-x1)/(2*h),(y2-y1)/(2*h)))
# (4) xdot = v/4 となる時刻
tt=sp.symbols('tt',positive=True)
X=sp.sqrt(16*lam**2+4*lam*v*tt)
print('(4) t =',sp.solve(sp.Eq(sp.diff(X,tt),v/4),tt))
# (5) 頂点 (0,-λ) での誕生：t=τ 後の x
tau=sp.symbols('tau',positive=True)
yy=-lam+v*tau; rr=lam+v*tau
print('(5) x =',sp.simplify(sp.sqrt(rr**2-yy**2)))
# (2) 負の y 軸と x 軸上の同位相点・逆位相点を波形で数値確認
T=1.0
def z(xx,yy,tt): 
    rr=np.hypot(xx,yy); return np.cos(2*np.pi*(tt/T-rr/L))+np.cos(2*np.pi*(tt/T-yy/L))
ts=np.linspace(0,1,201)
for yy in [-0.5,-1.0,-0.25,-0.75]:
    print(' y=',yy,' max|z| =',np.abs([z(0,yy,tt) for tt in ts]).max())
for xx in [1.0,2.0,0.5]:
    print(' x=',xx,' max|z| =',np.abs([z(xx,0,tt) for tt in ts]).max())
