# PHY-MOM-042 検証：2つの台をつなぐ糸の中央のおもり（中点の束縛＋水平運動量保存）
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
M,g,s,L=sp.symbols('M g s L',positive=True)
t=sp.symbols('t')
x1=sp.Function('x1')(t); th=sp.Function('th')(t)
x2=x1+2*s*sp.cos(th)
X=(x1+x2)/2; Y=L-s*sp.sin(th)
KE=sp.Rational(1,2)*M*x1.diff(t)**2+sp.Rational(1,2)*2*M*x2.diff(t)**2+sp.Rational(1,2)*M*(X.diff(t)**2+Y.diff(t)**2)
PE=M*g*Y
Lag=KE-PE
eqs=[sp.diff(Lag.diff(q.diff(t)),t)-Lag.diff(q) for q in (x1,th)]
a1,a2=sp.symbols('a1 a2'); v1,v2,q1,q2=sp.symbols('v1 v2 q1 q2')
rep={x1.diff(t,2):a1,th.diff(t,2):a2}
rep2={x1.diff(t):v1,th.diff(t):v2,x1:q1,th:q2}
E=[e.subs(rep).subs(rep2) for e in eqs]
sol=sp.solve(E,[a1,a2],dict=True)[0]
f=sp.lambdify((q1,q2,v1,v2,M,g,s,L),[sol[a1],sol[a2]])
Mv,gv,sv,Lv=1.0,9.8,1.0,2.0
def rhs(tt,y):
    A=f(y[0],y[1],y[2],y[3],Mv,gv,sv,Lv); return [y[2],y[3],A[0],A[1]]
ev=lambda tt,y: y[1]-np.pi/6; ev.terminal=False
ev2=lambda tt,y: y[1]-(np.pi/2-1e-6); ev2.terminal=True
r=solve_ivp(rhs,[0,5],[0,1e-7,0,1e-7*0],events=[ev,ev2],rtol=1e-11,atol=1e-13,dense_output=True)
# 初期条件: θ=0 から静かに（θ=0 は不安定なつり合いではなく重力で下がる：θ'' > 0 を確認）
print('initial accel',f(0,1e-7,0,0,Mv,gv,sv,Lv))
y6=r.y_events[0][0]
VP=y6[2]; VQ=VP+(-2*sv*np.sin(y6[1])*y6[3])
print('theta=pi/6: VP=',VP,' analytic',5/2*np.sqrt(gv*sv/23),'  VQ/VP=',VQ/VP)
# 張力：M x1''=T_P cosθ, 2M x2''=-T_Q cosθ
A=f(*y6,Mv,gv,sv,Lv)
x1dd=A[0]; th_,thd_=y6[1],y6[3]
x2dd=x1dd-2*sv*(np.cos(th_)*thd_**2+np.sin(th_)*A[1])
TP=Mv*x1dd/np.cos(th_); TQ=-2*Mv*x2dd/np.cos(th_)
print('T_P=',TP,' analytic',550*Mv*gv/529,'  T_Q/T_P=',TQ/TP)
yf=r.y_events[1][0]
print('near pi/2: VP=',yf[2],' analytic',5*np.sqrt(gv*sv/22))
X0=(0+0+2*sv)/2; Xf=(yf[0]+yf[0]+2*sv*np.cos(yf[1]))/2
print('R horizontal displacement',Xf-X0,' analytic s/4=',sv/4)
# 張力が常に正
ts=np.linspace(r.t[0]+1e-3,r.t[-1]-1e-3,400); mn=1e9
for tt in ts:
    y=r.sol(tt); A=f(*y,Mv,gv,sv,Lv); mn=min(mn,Mv*A[0]/np.cos(y[1]))
print('min T_P over motion',mn)
# エネルギー保存
def energy(y):
    x1v,thv,v1v,w=y; x2d=v1v-2*sv*np.sin(thv)*w; Xd=(v1v+x2d)/2; Yd=-sv*np.cos(thv)*w
    return 0.5*Mv*v1v**2+0.5*2*Mv*x2d**2+0.5*Mv*(Xd**2+Yd**2)+Mv*gv*(Lv-sv*np.sin(thv))
print('energy drift',energy(yf)-energy(r.y[:,0]))
