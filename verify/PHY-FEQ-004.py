# J011-1 検証：2つの粗い斜面にまたがるロープのつり合いの範囲・頂点の張力の不定性・滑り出した後の速さ
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
x,L,M,g=sp.symbols('x L M g',positive=True)
sa,ca,sb,cb,mu=sp.Rational(4,5),sp.Rational(3,5),sp.Rational(3,5),sp.Rational(4,5),sp.Rational(1,4)
lam=M/L
# 左（傾斜角 alpha）長さ x，右（beta）長さ L-x
TlL,TlU=lam*x*g*(sa-mu*ca), lam*x*g*(sa+mu*ca)
TrL,TrU=lam*(L-x)*g*(sb-mu*cb), lam*(L-x)*g*(sb+mu*cb)
xmin=sp.solve(sp.Eq(TlU,TrL),x)[0]; xmax=sp.solve(sp.Eq(TrU,TlL),x)[0]
print('xmin=',xmin,' xmax=',xmax, sp.N(xmin/L), sp.N(xmax/L))
# 一般式
al,be,m_=sp.symbols('alpha beta mu',positive=True)
xmin_g=sp.solve(sp.Eq(x*(sp.sin(al)+m_*sp.cos(al)),(L-x)*(sp.sin(be)-m_*sp.cos(be))),x)[0]
xmax_g=sp.solve(sp.Eq((L-x)*(sp.sin(be)+m_*sp.cos(be)),x*(sp.sin(al)-m_*sp.cos(al))),x)[0]
print('general xmin',sp.simplify(xmin_g),' xmax',sp.simplify(xmax_g))
# x=L/2 の張力の範囲
h={x:L/2}
lo=sp.Max(TlL.subs(h),TrL.subs(h)); hi=sp.Min(TlU.subs(h),TrU.subs(h)); print('T range at L/2:',sp.simplify(lo/(M*g)),sp.simplify(hi/(M*g)))
# 滑り出し（左へ）後: 動摩擦係数 mu で，M a = lam g [x(sa - mu ca) - (L-x)(sb + mu cb)]
F=lam*g*(x*(sa-mu*ca)-(L-x)*(sb+mu*cb))
W=sp.integrate(F,(x,xmax,L)); v2=sp.simplify(2*W/M); print('v^2 at B reaching apex =',v2,'=',sp.N(v2/(g*L)),'gL')
# 数値積分（xmax からわずかにずらして出発）
gn,Ln=9.8,1.0
def rhs(t,y):
    xx,vv=y; return [vv, gn/Ln*(xx*(0.8-0.25*0.6)-(Ln-xx)*(0.6+0.25*0.8))]
x0=float(xmax.subs(L,1))+1e-7
ev=lambda t,y: y[0]-Ln; ev.terminal=True
s=solve_ivp(rhs,[0,50],[x0,0],events=ev,rtol=1e-11,atol=1e-13)
print('numeric v^2/(gL)=',s.y_events[0][0][1]**2/(gn*Ln),' time',s.t_events[0][0])
