# O20081-1 : ばねにつながれた斜面台の上をすべり降りる小物体
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,k,g,th,t=sp.symbols('m M k g theta t',positive=True)
X,A,N=sp.symbols('X A N')
s,c=sp.sin(th),sp.cos(th)
# 台の加速度 A(右向き). 小物体の法線方向加速度 = A·n, n=(-s,c)
eqN=sp.Eq(m*(-A*s),N-m*g*c)
eqW=sp.Eq(M*A,-k*X+N*s)
sol=sp.solve([eqN,eqW],[A,N]); print('A=',sp.simplify(sol[A]),' N=',sp.simplify(sol[N]))
Om2=k/(M+m*s**2); Xc=m*g*s*c/k
print('check A=-Om2(X-Xc):',sp.simplify(sol[A]+Om2*(X-Xc)))
Om=sp.sqrt(Om2)
ell=sp.Rational(1,2)*g*s*sp.pi**2/Om2+2*Xc*c
print('ell=',sp.simplify(ell))
# 速さとエネルギー
v=g*s*sp.pi/Om
print('energy check:',sp.simplify(m*g*ell*s-(sp.Rational(1,2)*m*v**2+sp.Rational(1,2)*k*(2*Xc)**2)))
ratio=sp.simplify((sp.Rational(1,2)*k*(2*Xc)**2)/(sp.Rational(1,2)*m*v**2)); print('spring/KE=',ratio)
Nmin=sp.simplify(sol[N].subs(X,0)); Nend=sp.simplify(sol[N].subs(X,2*Xc)); print('N(0)=',sp.factor(Nmin),' N(end)=',sp.factor(Nend))
# 数値: 2自由度 (X, 斜面に沿った相対位置 sr) をラグランジュ方程式で積分
mv,Mv,kv,gv,tv=1.0,2.0,30.0,9.8,np.radians(35)
sv,cv=np.sin(tv),np.cos(tv)
ell_n=float(ell.subs({m:mv,M:Mv,k:kv,g:gv,th:tv})); Om_n=float(Om.subs({m:mv,M:Mv,k:kv,g:gv,th:tv}))
# 小物体位置 = (X + sr*cv, sr*sv) ; sr: 斜面下端からの距離
# L = 1/2 M Xd^2 + 1/2 m[(Xd+srd cv)^2+(srd sv)^2] - m g sr sv - 1/2 k X^2
def f(tt,y):
    Xv,sr,Xd,srd=y
    Mat=np.array([[Mv+mv, mv*cv],[mv*cv, mv]])
    rhs=np.array([-kv*Xv, -mv*gv*sv])
    acc=np.linalg.solve(Mat,rhs)
    return [Xd,srd,acc[0],acc[1]]
def ev(tt,y): return y[1]
ev.terminal=True
sl=solve_ivp(f,(0,10),[0,ell_n,0,0],events=ev,rtol=1e-11,atol=1e-12)
te=sl.t_events[0][0]; ye=sl.y_events[0][0]
print('reach foot t=',te,' pi/Om=',np.pi/Om_n,' X=',ye[0],' 2Xc=',2*float(Xc.subs({m:mv,k:kv,g:gv,th:tv})),' Xdot=',ye[2])
vb=np.hypot(ye[2]+ye[3]*cv, ye[3]*sv); print('block speed',vb,' theory',float(v.subs({m:mv,M:Mv,k:kv,g:gv,th:tv})))
