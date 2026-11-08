# O20141-1 検証：動ける台の円弧面から飛び出す小物体（傾いた垂直抗力と台の加速度）
from sympy import *
import numpy as np
from scipy.integrate import solve_ivp
M,m,R,g,a=symbols('M m R g alpha',positive=True)
th=symbols('theta')
# 理論式
u2=2*g*R*cos(a)*(M+m)/(M+m*sin(a)**2)
VD=-m*sqrt(u2)*cos(a)/(M+m)
ND=M*m*(u2/R+g*cos(a))/(M+m*sin(a)**2)
AD=ND*sin(a)/M
tanphi=(M+m)/M*tan(a)
s=u2*sin(2*a)/g
DX=-m/(M+m)*(R*(1+sin(a))+s)
print('N_D simplified:',factor(simplify(ND)))
print('s simplified:',simplify(s))
# 運動量・エネルギーの記号チェック（θ一般）
u=symbols('u',positive=True); V=symbols('V')
Vsol=solve(Eq(m*(V+u*cos(th))+M*V,0),V)[0]
KE=Rational(1,2)*m*((Vsol+u*cos(th))**2+(u*sin(th))**2)+Rational(1,2)*M*Vsol**2
usol=solve(Eq(KE,m*g*R*cos(th)),u)
assert any(simplify(us**2-u2.subs(a,th))==0 for us in usol)
# 地上から見た飛び出し角
vx=VD+sqrt(u2)*cos(a); vy=sqrt(u2)*sin(a)
assert simplify(vy/vx-tanphi)==0
# ラグランジュ方程式で数値積分（台の位置X，角θ：最下点から右回りに正）
Mn,mn,Rn,gn,an=2.0,1.0,1.0,9.8,np.pi/3
def rhs(t,y):
    X,T,Xd,Td=y
    # 物体位置 x=X+R sinT, y=-R cosT
    # 質量行列
    Mm=np.array([[Mn+mn, mn*Rn*np.cos(T)],[mn*Rn*np.cos(T), mn*Rn**2]])
    F=np.array([mn*Rn*np.sin(T)*Td**2, -mn*gn*Rn*np.sin(T)])
    acc=np.linalg.solve(Mm,F)
    return [Xd,Td,acc[0],acc[1]]
ev=lambda t,y: y[1]-an; ev.terminal=True; ev.direction=1
sol=solve_ivp(rhs,[0,20],[0,-np.pi/2,0,0],events=ev,rtol=1e-11,atol=1e-12,max_step=1e-3)
X,T,Xd,Td=sol.y_events[0][0]
acc=rhs(0,[X,T,Xd,Td])
Xdd,Tdd=acc[2],acc[3]
# 物体の加速度
ax=Xdd+Rn*(np.cos(T)*Tdd-np.sin(T)*Td**2); ay=Rn*(np.sin(T)*Tdd+np.cos(T)*Td**2)
# 垂直抗力：m a = N n + m g  (n=内向き(-sin,cos))
n=np.array([-np.sin(T),np.cos(T)])
Nnum=np.dot(mn*np.array([ax,ay+gn]),n)
sub={M:Mn,m:mn,R:Rn,g:gn,a:an}
print('u_D num',Rn*Td,'th',float(sqrt(u2.subs(sub))))
print('V_D num',Xd,'th',float(VD.subs(sub)))
print('N_D num',Nnum,'th',float(ND.subs(sub)))
print('A_D num',Xdd,'th',float(AD.subs(sub)))
print('X at D num',X,'th',-mn/(Mn+mn)*Rn*(1+np.sin(an)))
assert abs(Rn*Td-float(sqrt(u2.subs(sub))))<1e-6
assert abs(Xd-float(VD.subs(sub)))<1e-6
assert abs(Nnum-float(ND.subs(sub)))<1e-5
assert abs(Xdd-float(AD.subs(sub)))<1e-5
# 接線方向の力のつり合いチェック（垂直抗力以外の力はない）
tang=np.array([np.cos(T),np.sin(T)])
assert abs(np.dot(mn*np.array([ax,ay+gn]),tang))<1e-6
# 飛行：台は等速，物体は放物運動．DE面（D と同じ高さ）に戻る時刻
uD=Rn*Td; tf=2*uD*np.sin(an)/gn
xo=X+Rn*np.sin(an)+(Xd+uD*np.cos(an))*tf
Xb=X+Xd*tf
sn=xo-(Xb+Rn*np.sin(an))
print('s num',sn,'th',float(s.subs(sub)))
print('DX num',Xb,'th',float(DX.subs(sub)))
assert abs(sn-float(s.subs(sub)))<1e-6 and abs(Xb-float(DX.subs(sub)))<1e-6
# 弧の途中で N>0 か確認（数値解の各点で）
Nmin=1e9
for i in range(len(sol.t)):
    y=sol.y[:,i]; d=rhs(0,y); T_,Td_=y[1],y[3]
    ax_=d[2]+Rn*(np.cos(T_)*d[3]-np.sin(T_)*Td_**2); ay_=Rn*(np.sin(T_)*d[3]+np.cos(T_)*Td_**2)
    Nmin=min(Nmin,np.dot(mn*np.array([ax_,ay_+gn]),[-np.sin(T_),np.cos(T_)]))
print('min N along arc',Nmin); assert Nmin>=-1e-9
print('OK')
