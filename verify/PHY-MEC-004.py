# Q105-1 軽い滑車にかけた鎖（端におもり）が滑車から離れる瞬間
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
M,m,L,g,x=sp.symbols('M m L g x',positive=True)
lam=M/L
x0=sp.solve(sp.Eq(lam*x+m,lam*(L-x)),x)[0]; print('x0',sp.simplify(x0))
a=g*(lam*(2*x-L)+m)/(M+m)
# エネルギー保存 (x0 から静かに放す): v^2
U=-(lam*x**2/2+lam*(L-x)**2/2)*g - m*g*x
v2=sp.simplify(2*(U.subs(x,x0)-U)/(M+m)); print('v^2',sp.factor(v2))
T=sp.simplify(lam*(L-x)*(g+a)); TA=sp.simplify((lam*x+m)*(g-a)); print('T consistency',sp.simplify(T-TA))
print('T',sp.factor(T))
sols=sp.solve(sp.Eq(T,lam*v2),x); print('x_sep',sols)
for xs in sols:
    print(' v_sep^2',sp.factor(sp.simplify(v2.subs(x,xs))), ' L-x',sp.simplify(L-xs))
# 数値積分: 少しずらして放し，T=lam v^2 の瞬間を探す
for (Mv,mv,Lv) in [(2.0,0.0,1.0),(2.0,0.6,1.0),(3.0,2.5,2.0)]:
    gv=9.8; lv=Mv/Lv; xe=(Lv-mv/lv)/2
    def rhs(t,y):
        X,V=y; return [V, gv*(lv*(2*X-Lv)+mv)/(Mv+mv)]
    def ev(t,y):
        X,V=y; A=gv*(lv*(2*X-Lv)+mv)/(Mv+mv); return lv*(Lv-X)*(gv+A)-lv*V**2
    ev.terminal=True
    so=solve_ivp(rhs,[0,50],[xe+1e-7,0],events=ev,rtol=1e-12,atol=1e-14)
    vs=so.y_events[0][0][1]; xs=so.y_events[0][0][0]
    print('M,m,L',Mv,mv,Lv,' v sim',vs,' theory',0.5*np.sqrt(gv*(Mv+mv)*Lv/Mv),' x sim',xs,' theory',xe+(Mv+mv)*Lv/(2*np.sqrt(2)*Mv))
