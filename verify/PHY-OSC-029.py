import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
g=r=1.0
# (2) symbolic cos thF = 2/3 cos thV
v0,th,thV=sp.symbols('v0 theta thetaV')
v2=v0**2-2*g*r*(1-sp.cos(th))
F=v2/r+g*sp.cos(th)   # per unit mass, positive outward rail
sol=sp.solve(sp.Eq(F.subs(v0**2,2*g*r*(1-sp.cos(thV))),0),sp.cos(th))
print('cos thF =',sol, ' expected 2/3 cos thV')
# (3) leave at cos thF=-c, pass through centre
def miss(c):
    thF=np.pi-np.arccos(c); v=np.sqrt(g*r*c)
    x0,y0=r*np.sin(thF),-r*np.cos(thF); vx,vy=v*np.cos(thF),v*np.sin(thF)
    t=x0/(-vx)  # time to reach x=0
    return y0+vy*t-0.5*g*t*t
c=brentq(miss,0.1,0.99); print('c=',c,'1/sqrt3=',1/np.sqrt(3))
cosV=-1.5*c; print('thetaV deg',np.degrees(np.arccos(cosV)),' v0^2=',2*g*r*(1-cosV),' 2+sqrt3=',2+np.sqrt(3))
thF=np.pi-np.arccos(c); v=np.sqrt(g*r*c)
x0,y0=r*np.sin(thF),-r*np.cos(thF); vx,vy=v*np.cos(thF),v*np.sin(thF)
hmax=y0+vy**2/(2*g)+r; print('hmax above A',hmax,' 1+4/(3sqrt3)=',1+4/(3*np.sqrt(3)),' h2=',1+np.sqrt(3)/2)
# where it hits circle again
f=lambda t:(x0+vx*t)**2+(y0+vy*t-0.5*g*t*t)**2-r**2
t2=brentq(f,0.5,5); X,Y=x0+vx*t2,y0+vy*t2-0.5*g*t2*t2
print('hit point',X,Y,'angle from bottom (deg, ccw):',np.degrees(np.arctan2(X,-Y))%360)
# numerical check: full motion in outer rail using ODE with constraint until F=0
def rhs(t,s):
    th,om=s; return [om,-g/r*np.sin(th)]
ev=lambda t,s: (r*s[1]**2+g*np.cos(s[0]))  # F/m
ev.terminal=True; ev.direction=-1
s=solve_ivp(rhs,[0,20],[0,np.sqrt(2+np.sqrt(3))],events=ev,rtol=1e-11,atol=1e-12)
print('ODE leave angle',np.degrees(s.y_events[0][0][0]),' expected',np.degrees(thF))
# landing on the floor level (y=-r), left of A
fl=lambda t: y0+vy*t-0.5*g*t*t + r
tl=brentq(fl,0.5,5); print('land x =',x0+vx*tl)
cs=sp.Rational(1,1)/sp.sqrt(3); thFs=sp.pi-sp.acos(cs); vs=sp.sqrt(cs)
T=sp.symbols('T',positive=True)
X0,Y0=sp.sin(thFs),-sp.cos(thFs); VX,VY=vs*sp.cos(thFs),vs*sp.sin(thFs)
sols=sp.solve(sp.Eq(Y0+VY*T-T**2/2,-1),T)
for s_ in sols:
    if s_.is_positive: print('exact land x =',sp.nsimplify(sp.simplify(X0+VX*s_)), sp.N(X0+VX*s_))
