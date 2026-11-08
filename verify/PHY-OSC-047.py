# Q082-1: 放物線形の丘をすべる物体：放物運動との比較で N・R = 一定
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
x,h,d,g,m,w=sp.symbols('x h d g m w',positive=True)
y=h*(1-4*x**2/d**2)
yp=sp.diff(y,x); ypp=sp.diff(y,x,2)
R=sp.simplify((1+yp**2)**sp.Rational(3,2)/sp.Abs(ypp))
u2=g*d**2/(8*h)
cos=1/sp.sqrt(1+yp**2)
assert sp.simplify(R-u2/(g*cos**3))==0
# 水平投射が放物線に一致：y = h - g x^2/(2u^2)
assert sp.simplify(h-g*x**2/(2*u2)-y)==0
v2=w**2+2*g*(h-y)
N=m*g*cos-m*v2/R
print('N*R =',sp.simplify(N*R)); assert sp.simplify(N*R-m*(u2-w**2))==0
vals={h:9.0,d:60.0,g:9.8,m:50.0,w:10.0}
print('u^2=',float(u2.subs(vals)),' u=',float(sp.sqrt(u2).subs(vals)))
for xv in [0,15,30]:
    print('x=',xv,' R=',float(R.subs(vals).subs(x,xv)),' N=',float(N.subs(vals).subs(x,xv)))
# 数値：放物線に拘束した運動を直接積分し，拘束力を計算（ラグランジュ形式で x の運動方程式）
hv,dv,gv,mv,wv=9.0,60.0,9.8,50.0,10.0
f=lambda X: hv*(1-4*X**2/dv**2); fp=lambda X:-8*hv*X/dv**2; fpp=-8*hv/dv**2
def rhs(t,s):
    X,Xd=s
    # (1+f'^2) Xdd + f' f'' Xd^2 = -g f'
    Xdd=(-gv*fp(X)-fp(X)*fpp*Xd**2)/(1+fp(X)**2)
    return [Xd,Xdd]
ev=lambda t,s: s[0]-30; ev.terminal=True
sol=solve_ivp(rhs,[0,20],[0,wv],events=ev,rtol=1e-12,atol=1e-12,dense_output=True)
for X_target in [0,15,29.999]:
    # 時刻を探す
    ts=np.linspace(0,sol.t[-1],200001); Xs=sol.sol(ts)[0]; i=np.argmin(abs(Xs-X_target))
    X,Xd=sol.sol(ts[i]); Xdd=rhs(0,[X,Xd])[1]
    ax=Xdd; ay=fpp*Xd**2+fp(X)*Xdd
    Fx=mv*ax; Fy=mv*ay+mv*gv     # 垂直抗力（ベクトル）
    print('num N at x=%.2f: %.4f'%(X,np.hypot(Fx,Fy)))
print('OK')
