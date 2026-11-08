# PHY-ENG-014 検証: 一定の力で引かれる A と糸でつながる B,C（V字）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
t=sp.symbols('t'); M,m,L,F=sp.symbols('M m L F',positive=True)
x=sp.Function('x')(t); p=sp.Function('p')(t)
xd,pd=x.diff(t),p.diff(t)
# B=(x-L cos p, L sin p), C=(x-L cos p,-L sin p)
T=sp.Rational(1,2)*M*xd**2+m*((xd+L*sp.sin(p)*pd)**2+(L*sp.cos(p)*pd)**2)
V=-F*x
Lg=T-V
eqs=[sp.diff(Lg.diff(q.diff(t)),t)-Lg.diff(q) for q in (x,p)]
a1,a2=sp.symbols('a1 a2')
sol=sp.solve([e.subs({x.diff(t,2):a1,p.diff(t,2):a2}) for e in eqs],[a1,a2],dict=True)[0]
X,P,XD,PD=sp.symbols('X P XD PD')
rep={xd:XD,pd:PD,x:X,p:P}
f1=sp.lambdify((X,P,XD,PD,M,m,L,F),sol[a1].subs(rep))
f2=sp.lambdify((X,P,XD,PD,M,m,L,F),sol[a2].subs(rep))
def run(Mv,mv,Lv,Fv,th0):
    def rhs(tt,y): return [y[2],y[3],f1(*y,Mv,mv,Lv,Fv),f2(*y,Mv,mv,Lv,Fv)]
    ev=lambda tt,y:y[1]; ev.terminal=True; ev.direction=-1
    s=solve_ivp(rhs,[0,100],[0,th0,0,0],events=ev,rtol=1e-11,atol=1e-12,dense_output=True)
    ts=s.t_events[0][0]; yE=s.y_events[0][0]
    # 張力: A の運動方程式 M a_A = F - 2 T cos p
    Ts=[]
    for tt in np.linspace(0,ts,400):
        y=s.sol(tt); aA=f1(*y,Mv,mv,Lv,Fv); Ts.append((Fv-Mv*aA)/(2*np.cos(y[1])))
    return ts,yE,min(Ts),Ts[-1],s
for (Mv,mv,Lv,Fv,th0) in [(2,1,1,1,np.pi/3),(1,1,1,3,1.2),(3,0.5,2,1,0.5)]:
    ts,yE,Tmin,Tend,s=run(Mv,mv,Lv,Fv,th0)
    w=abs(Lv*yE[3]); u=yE[2]
    w2=2*Fv*Lv*(1-np.cos(th0))/(Mv+2*mv)
    xA=Fv*ts**2/(2*(Mv+2*mv))+2*mv*Lv*(1-np.cos(th0))/(Mv+2*mv)
    Tf=mv*(Mv*w2/Lv+Fv)/(Mv+2*mv)
    print('w^2 num/formula',w**2,w2,' u num/formula',u,Fv*ts/(Mv+2*mv),' xA',yE[0],xA,' T_end',Tend,Tf,' Tmin',Tmin)
    # 弾性衝突後: p' -> -p'、最大半角が th0 に戻るか
    y0=yE.copy(); y0[3]=-y0[3]; y0[1]=1e-12
    def rhs(tt,y): return [y[2],y[3],f1(*y,Mv,mv,Lv,Fv),f2(*y,Mv,mv,Lv,Fv)]
    ev=lambda tt,y:y[3]; ev.terminal=True; ev.direction=-1
    s2=solve_ivp(rhs,[0,100],y0,events=ev,rtol=1e-11,atol=1e-12)
    print('  max half-angle after elastic',s2.y_events[0][0][1],'th0',th0)
# 記号計算：重心系の位置エネルギーで w^2 を別導出
th0=sp.symbols('theta0',positive=True); w=sp.symbols('w',positive=True)
# 相対運動エネルギー: A の重心に対する x 速度と B,C の y 速度（衝突時 B,C の x 速度=A と等しい → 相対KE = m w^2）
print(sp.solve(sp.Eq(m*w**2, F*2*m*L*(1-sp.cos(th0))/(M+2*m)),w))
