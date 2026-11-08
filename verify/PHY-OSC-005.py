# M30-1: smooth insulating incline angle a; in-plane coords x (horizontal), y (up-slope). Uniform E along +x with qE = m g sin a.
# ball charge q>0 on string length r from O, circles counterclockwise (seen from above the plane). theta from +x axis.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,g,a,r,v0,th=sp.symbols('m g alpha r v_0 theta',positive=True)
qE=m*g*sp.sin(a)
U=m*g*sp.sin(a)*r*sp.sin(th)-qE*r*sp.cos(th)     # potential energy (gravity along slope + field)
v2=v0**2-2*(U-U.subs(th,0))/m
# radial: T + (radial component of in-plane force toward outside?) in-plane force F=(qE, -mg sin a); outward radial unit (cos th, sin th)
Frad=qE*sp.cos(th)-m*g*sp.sin(a)*sp.sin(th)
T=sp.simplify(m*v2/r + Frad)   # inward net force T - Frad = m v^2/r
print('T(theta)=',T)
print('dT/dth at 135deg =',sp.simplify(sp.diff(T,th).subs(th,3*sp.pi/4)), ' T(135)-T(th) >=0 check below')
Tmin=sp.simplify(T.subs(th,3*sp.pi/4)); print('T at 135deg =',Tmin)
assert sp.simplify(Tmin-(m*v0**2/r-m*g*sp.sin(a)*(3*sp.sqrt(2)+2)))==0
v0min=sp.solve(sp.Eq(Tmin,0),v0)[0]; print('v0min=',v0min)
# Method 2: effective gravity g' = sqrt2 g sin a; effective top at 135deg; v_top^2 >= g' r ; height diff top vs theta=0 in effective sense: r + r cos45 -> drop
gp=sp.sqrt(2)*g*sp.sin(a)
v0min2=sp.sqrt(gp*r+2*gp*(r+r/sp.sqrt(2)))
assert sp.simplify(v0min**2-v0min2**2)==0
# numeric: integrate circular motion with T, then cut at 135deg with min speed
gv,av,rv=9.8,np.radians(30),1.0; gs=gv*np.sin(av)
v0v=float(v0min.subs({g:gv,a:av,r:rv}))
def f(t,y):  # y=[theta, thetadot]
    th_,w_=y
    Ft=-gs*np.sin(th_)*1 - gs*np.cos(th_)*1  # tangential accel: F.(-sin,cos)/m = gs*(-sin) + (-gs)*cos
    return [w_, Ft/rv]
def ev(t,y): return y[0]-3*np.pi/4
ev.terminal=True
s=solve_ivp(f,[0,5],[0,v0v/rv],events=ev,rtol=1e-12,atol=1e-12)
th_,w_=s.y_events[0][0]
Tm=rv*w_**2+(gs*np.cos(th_)-gs*np.sin(th_)); print('T/m at 135 (numeric)=',Tm)
assert abs(Tm)<1e-6
# projectile after cut
p=rv*np.array([np.cos(th_),np.sin(th_)]); vv=rv*w_*np.array([-np.sin(th_),np.cos(th_)])
def fp(t,y): return [y[2],y[3],gs,-gs]
ev2=lambda t,y: y[1]-y[0]; ev2.terminal=True; ev2.direction=-1
s2=solve_ivp(fp,[0,5],[*p,*vv],events=ev2,rtol=1e-12,atol=1e-12)
P=s2.y_events[0][0][:2]; print('cross y=x at',P,' dist=',np.hypot(*P),' sqrt2 r=',np.sqrt(2)*rv)
assert abs(np.hypot(*P)-np.sqrt(2)*rv)<1e-8
# Also check: ball passes through effective lowest? and that after cut ball stays on incline (N=mg cos a >0) trivially
print('OK')
# check T minimum at 135deg over theta grid
Tf=sp.lambdify(th,(T/m).subs({g:9.8,a:0.5,r:1.0,v0:5.0}))
grid=np.linspace(0,2*np.pi,100001); i=np.argmin(Tf(grid)); print('argmin theta(deg)=',np.degrees(grid[i]))
assert abs(np.degrees(grid[i])-135)<0.01
print('OK2')
# speed at S vs sqrt(3*sqrt2*g*r*sin a)
vS=np.hypot(*s2.y_events[0][0][2:]); print('speed at S=',vS,' formula=',np.sqrt(3*np.sqrt(2)*gv*rv*np.sin(av)))
assert abs(vS-np.sqrt(3*np.sqrt(2)*gv*rv*np.sin(av)))<1e-8
# alternative form of T
assert sp.simplify(T-(m*v0**2/r+m*g*sp.sin(a)*(3*sp.cos(th)-3*sp.sin(th)-2)))==0
print('OK3')
