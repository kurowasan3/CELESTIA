# P18-1: 近日点・遠日点の曲率半径（調和平均）と近日点での衝突後の軌道
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
G,Ms,r1,r2,v1,v2,W,M,m,u=sp.symbols('G M_s r_1 r_2 v_1 v_2 W M m u',positive=True)
s=sp.solve([sp.Eq(r1*v1,r2*v2), sp.Eq(v1**2/2-G*Ms/r1, v2**2/2-G*Ms/r2)],[v1,v2],dict=True)[0]
rho1=sp.simplify((v1**2*r1**2/(G*Ms)).subs(s)); rho2=sp.simplify((v2**2*r2**2/(G*Ms)).subs(s))
print('rho1',rho1,' rho2',rho2, ' equal:',sp.simplify(rho1-rho2)==0)
print('1/rho - (1/r1+1/r2)/2 =',sp.simplify(1/rho1-(1/r1+1/r2)/2))
# 衝突後: rho' = W^2 r1^2/(G Ms), もう一方の頂点 r'
rhop=W**2*r1**2/(G*Ms)
rp=sp.simplify(1/(2/rhop-1/r1)); print("r' =",rp)
# エネルギー・面積速度で独立に解いて照合
x=sp.symbols('x',positive=True); Wx=sp.symbols('Wx',positive=True)
sol=sp.solve([sp.Eq(r1*W,x*Wx), sp.Eq(W**2/2-G*Ms/r1,Wx**2/2-G*Ms/x)],[x,Wx],dict=True)
print('energy/area solutions:',[sp.simplify(q[x]) for q in sol])
# 円軌道になる u
Wc=sp.sqrt(G*Ms/r1); uc=sp.solve(sp.Eq((M*v1-m*u)/(M+m),Wc),u)[0]; print('u circ',uc)
# 数値: 軌道を積分して曲率半径と r' を確認
GM=1.0; R1,R2=1.0,2.0; V1=np.sqrt(2*GM*R2/(R1*(R1+R2)))
def f(t,y):
    r=np.hypot(y[0],y[1]); return [y[2],y[3],-GM*y[0]/r**3,-GM*y[1]/r**3]
def curv(y):
    a=f(0,y); vx,vy=y[2],y[3]; ax,ay=a[2],a[3]; return (vx*vx+vy*vy)**1.5/abs(vx*ay-vy*ax)
y0=[R1,0,0,V1]; ev=lambda t,y:y[1]; ev.direction=-1; ev.terminal=True
r=solve_ivp(f,[1e-6,100],y0,events=ev,rtol=1e-12,atol=1e-12)
print('rho at peri',curv(np.array(y0)),'rho at aph',curv(r.y[:,-1]),'harmonic',2*R1*R2/(R1+R2),'aph dist',-r.y[0,-1])
# 衝突で速さを 0.8 V1 に
Wn=0.8*V1; y0=[R1,0,0,Wn]; r=solve_ivp(f,[1e-6,100],y0,events=ev,rtol=1e-12,atol=1e-12)
print("r' num",-r.y[0,-1],'formula',float(rp.subs({W:Wn,r1:R1,G:1,Ms:1})), ' rho\'',Wn**2, '<r1 ->A becomes aphelion')
