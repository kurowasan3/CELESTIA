# Q094-1 球のバックスピン・動ける板
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,r,v0,w0,mu,g=sp.symbols('m M r v0 omega0 mu g',positive=True)
I=sp.Rational(2,5)*m*r**2
s0=v0+w0*r
k=1/m+r**2/I+1/M
J=-s0/k                     # 摩擦の力積（右向き正）
t1=sp.simplify(-J/(mu*m*g))
v1=sp.simplify(v0+J/m); V1=sp.simplify(-J/M); w1=sp.simplify(-w0 - J*r/I)
print('t1',sp.factor(t1)); print('v1',sp.factor(v1)); print('V1',sp.factor(V1)); print('w1',sp.factor(w1))
# 転がり条件確認
print('roll check',sp.simplify(v1-w1*r-V1))
# 角運動量（床上の固定点）保存の別解
L0=m*r*v0+I*(-w0); L1=m*r*v1+I*w1
print('L cons',sp.simplify(L0-L1))
# 運動量保存
print('p cons',sp.simplify(m*v0-(m*v1+M*V1)))
Q=sp.simplify(sp.Rational(1,2)*m*v0**2+sp.Rational(1,2)*I*w0**2-(sp.Rational(1,2)*m*v1**2+sp.Rational(1,2)*I*w1**2+sp.Rational(1,2)*M*V1**2))
print('Q',sp.factor(Q), ' target', sp.simplify(Q-m*M*s0**2/(7*M+2*m)))
cond=sp.solve(sp.Eq(v1,0),w0); print('v1=0 at w0=',cond)
# 数値積分
vals=dict(m=1.0,M=3.0,r=0.1,v0=2.0,w0=60.0,mu=0.3,g=9.8)
def f(t,y):
    v,w,V=y; s=v-w*vals['r']-V
    fr=-vals['mu']*vals['m']*vals['g']*np.tanh(s/1e-6)
    Iv=0.4*vals['m']*vals['r']**2
    return [fr/vals['m'], -fr*vals['r']/Iv, -fr/vals['M']]
sol=solve_ivp(f,[0,3],[vals['v0'],-vals['w0'],0],max_step=1e-3,rtol=1e-9,atol=1e-12,method='LSODA')
sub={m:1,M:3,r:sp.Rational(1,10),v0:2,w0:60,mu:sp.Rational(3,10),g:sp.Rational(49,5)}
print('num v,w,V',sol.y[:,-1],' theory',[float(v1.subs(sub)),float(w1.subs(sub)),float(V1.subs(sub))])
print('t1 theory',float(t1.subs(sub)))
# (4) 任意の板の運動: 板速度 U(t) を与える
def f2(t,y):
    v,w=y; U=3*np.sin(2*t)*np.exp(-0.3*t) if t<8 else 0.0
    s=v-w*0.1-U; fr=-0.3*9.8*np.tanh(s/1e-6)
    return [fr, -fr*0.1/(0.4*0.01)]
s2=solve_ivp(f2,[0,15],[0,0],max_step=1e-3,method='LSODA')
print('(4) final v,w',s2.y[:,-1])
