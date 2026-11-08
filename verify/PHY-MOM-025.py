# P4 類題：折り返した鎖の片端を放す（折り返し点の速さ v/2，支点の力）
import sympy as sp, numpy as np
M,L,g,x,t=sp.symbols('M L g x t',positive=True)
lam=M/L
v=sp.sqrt(2*g*x)
lf=(L+x)/2; lm=(L-x)/2
p=lam*lm*v
# 時間微分： dx/dt = v
dpdt=sp.simplify(sp.diff(p,x)*v)
F=sp.simplify(M*g-dpdt)
print('F=',sp.factor(F))
assert sp.simplify(F-M*g*(L+3*x)/(2*L))==0
# 別解：静止部分の重さ + 運動量流入 (dm/dt)*v ; dm/dt = lam * v/2
F2=lam*g*lf+lam*(v/2)*v
assert sp.simplify(F-F2)==0
xs=sp.solve(sp.diff(p,x),x); print('p max at',xs, ' F there',sp.simplify(F.subs(x,xs[0])))
print('F(0)=',F.subs(x,0),' F(L)=',sp.simplify(F.subs(x,L)))
# エネルギー損失
U0=-M*g*L/4; Uf=-M*g*L/2
print('lost energy=',sp.simplify(U0-Uf))
tt=sp.symbols('tt',positive=True)
loss=sp.integrate(lam*(g*tt)**3/4,(tt,0,sp.sqrt(2*L/g)))
print('integral of loss rate=',sp.simplify(loss)); assert sp.simplify(loss-M*g*L/4)==0
# 力学的エネルギーの確認：x における KE+PE と初期との差 = 積分損失
KE=sp.Rational(1,2)*lam*lm*v**2
# PE: 静止部 長さ lf 重心深さ lf/2；動く部 上端深さ x 下端深さ lf 重心深さ (x+lf)/2
PE=-lam*g*(lf*lf/2+lm*(x+lf)/2)
E=sp.simplify(KE+PE); print('E(x)=',sp.expand(E))
lossx=sp.integrate(lam*(g*tt)**3/4,(tt,0,sp.sqrt(2*x/g)))
assert sp.simplify(E-(U0-lossx))==0
print('energy bookkeeping OK')
# 数値チェック：p(t) を差分して F を確認
Mn,Ln,gn=2.0,1.5,9.8
ts=np.linspace(0.01,np.sqrt(2*Ln/gn)*0.99,500)
xn=gn*ts**2/2; pn=(Mn/Ln)*(Ln-xn)/2*gn*ts
Fn=Mn*gn-np.gradient(pn,ts)
print('max dev',np.max(abs(Fn-Mn*gn*(Ln+3*xn)/(2*Ln))[2:-2]))
