# P24 類題：水漏れする缶の振り子（容器＋液体の重心が最も低いのは重心が液面と一致するとき）
import sympy as sp, numpy as np
M,b,rho,S,x=sp.symbols('M b rho S x',positive=True)
yG=(M*b+rho*S*x**2/2)/(M+rho*S*x)
dy=sp.simplify(sp.diff(yG,x)); print('dyG/dx numerator:',sp.factor(sp.numer(sp.together(dy))))
# 極小条件 rho S x^2/2 + M x - M b = 0 と yG = x の同値性
xs=sp.solve(sp.Eq(rho*S*x**2/2+M*x-M*b,0),x)[0]
print('x*=',xs,' yG(x*)-x* =',sp.simplify(yG.subs(x,xs)-xs))
print('yG(2b)-b =',sp.simplify(yG.subs(x,2*b)-b))
vals={M:400,b:10,rho:1,S:50}
xn=float(xs.subs(vals)); print('x* =',xn,'cm')
g=9.8; l=1.00
for xx in [20,xn,0]:
    y=float(yG.subs(vals).subs(x,xx)); print('x=',xx,' yG=',y,' T=',2*np.pi*np.sqrt((l-y/100)/g))
# 数値で最小を確認
xx=np.linspace(0,20,200001); yy=np.array([(4000+25*t**2)/(400+50*t) for t in xx])
print('numeric min at',xx[yy.argmin()],yy.min())
