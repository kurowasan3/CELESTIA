# M18-1: pendulum in a box on smooth floor released from horizontal, string cut at lowest point
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,l,H=sp.symbols('m M g l H',positive=True)
v,V=sp.symbols('v V',positive=True)
s=sp.solve([sp.Eq(m*v,M*V),sp.Eq(m*g*l,m*v**2/2+M*V**2/2)],[v,V],dict=True)[0]
print('(2)',sp.simplify(s[v]),sp.simplify(s[V]))
xrel=sp.simplify((s[v]+s[V])*sp.sqrt(2*H/g))
print('(3) x_rel',sp.simplify(xrel), ' compare', sp.simplify(xrel-2*sp.sqrt(l*H*(m+M)/M)))
# box position at landing by direct kinematics: D1 + V t
D1=m*l/(m+M)
Xland=sp.simplify(D1+s[V]*sp.sqrt(2*H/g))
print('(3) X_land direct',Xland,' CM form diff',sp.simplify(Xland-m*(l+xrel)/(m+M)))
Hs=sp.solve(sp.Eq(Xland-xrel,0),H)
print('(4) H=',Hs, ' expected', m**2*l/(4*M*(m+M)))
# numeric
mm,MM,gg,ll,HH=0.5,2.0,9.8,0.8,0.6
def rhs(t,y):
    X,th,Xd,thd=y
    A=np.array([[MM+mm, mm*ll*np.cos(th)],[mm*ll*np.cos(th), mm*ll**2]])
    b=np.array([mm*ll*np.sin(th)*thd**2, -mm*gg*ll*np.sin(th)])
    Xdd,thdd=np.linalg.solve(A,b); return [Xd,thd,Xdd,thdd]
ev=lambda t,y:y[1]; ev.terminal=True; ev.direction=-1
r=solve_ivp(rhs,[0,5],[0,np.pi/2,0,0],events=ev,rtol=1e-12,atol=1e-13,max_step=1e-4)
X,th,Xd,thd=r.y
# tension: along string toward hook, in box frame: T - m g cos th - m Xdd sin th ... compute from ball accel in ground frame
Xdd=np.array([rhs(0,r.y[:,i])[2] for i in range(len(r.t))]); thdd=np.array([rhs(0,r.y[:,i])[3] for i in range(len(r.t))])
ax=Xdd+ll*(np.cos(th)*thdd-np.sin(th)*thd**2); ay=ll*(np.sin(th)*thdd+np.cos(th)*thd**2)
# string unit from ball to hook: (-sin th, cos th)
T=mm*(ax*(-np.sin(th))+ (ay+gg)*np.cos(th))
print('min tension (after start)',T[5:].min(),' at bottom',T[-1])
print('D at bottom',X[-1],' th',mm*ll/(mm+MM))
vb=Xd[-1]+ll*thd[-1]; print('ball v',vb,' box V',Xd[-1],' th',-np.sqrt(2*MM*gg*ll/(mm+MM)),mm/MM*np.sqrt(2*MM*gg*ll/(mm+MM)))
t=np.sqrt(2*HH/gg); Xl=X[-1]+Xd[-1]*t; xb=vb*t   # ball ground x (starts at 0 at bottom since ball at hook x)
print('x_rel num',Xl-(X[-1]+xb),' th',2*np.sqrt(ll*HH*(mm+MM)/MM))
print('ball ground landing x',X[-1]+xb,' (CM check: m*xb+M*X should equal m*l):',mm*(X[-1]+xb)+MM*Xl, mm*ll)
# alternate: D - v t = 0
Hs2=sp.solve(sp.Eq(D1-s[v]*sp.sqrt(2*H/g),0),H); print('(4) alt H=',Hs2)
