# PHY-ESF-026 検証: 一様帯電球（トムソン型原子）の電子を一様電場で引き出す
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
k,e,R,E,x=sp.symbols('k e R E x',positive=True)
K=k*e**2/R**3
Uin=K*x**2/2; C=sp.symbols('C'); Cs=sp.solve(sp.Eq(Uin.subs(x,R),-k*e**2/R+C),C)[0]
Uout=-k*e**2/x+Cs; print('C=',Cs)
# 急に電場: 全エネルギー0, 最大変位 2eE/K
E1=sp.solve(sp.Eq(2*e*E/K,R),E)[0]; print('E1=',E1)
Ut=Uout-e*E*x; xs=sp.solve(sp.diff(Ut,x),x)[0]; Umax=sp.simplify(Ut.subs(x,xs))
E2=sp.solve(sp.Eq(Umax,0),E)[0]; print('x*=',xs,' Umax=',Umax,' E2=',E2,' x*(E2)/R=',sp.simplify(xs.subs(E,E2)/R))
E3=sp.solve(sp.Eq(e*E/K,R),E)[0]; print('E3 (slow)=',E3)
Ev=sp.Rational(27,50)*k*e/R**2
xm=sp.solve(sp.Eq(Ut.subs(E,Ev),0),x); print('turning points with E=27/50 ke/R^2:',[sp.simplify(r/R) for r in xm])
# 数値（k=e=m=R=1）
def force(xx,Ev):  # 電子に働く力（+x 向きを正）
    return (-xx if abs(xx)<=1 else -np.sign(xx)/xx**2)+Ev
def run(Ev,T=200):
    s=solve_ivp(lambda t,y:[y[1],force(y[0],Ev)],[0,T],[0,0],rtol=1e-10,atol=1e-12,max_step=0.01)
    return s.y[0].max()
for Ev in [0.49,0.51,0.54,0.56,0.5625-1e-3,0.5625+1e-3]:
    print('E=',Ev,' max x=',run(Ev))
# ゆっくり増やす: E=eps*t
for eps in [1e-3,3e-4]:
    def rhs(t,y): return [y[1],force(y[0],eps*t)]
    ev=lambda t,y:y[0]-3; ev.terminal=True
    s=solve_ivp(rhs,[0,2/eps],[0,0],events=ev,rtol=1e-10,atol=1e-12,max_step=0.05)
    # 球の外へ出た瞬間の電場
    xs_=s.y[0]; idx=np.argmax(xs_>1.0)
    print('slow ramp eps',eps,' E when leaving sphere',eps*s.t[idx],' E at escape(x=3)',eps*s.t_events[0][0])
