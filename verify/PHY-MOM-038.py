# S0110-1: launch on rough floor, non-bouncing landing with friction impulse, slide; maximize total distance
import sympy as sp, numpy as np
v,g,mu,th = sp.symbols('v g mu theta', positive=True)
u = v*(sp.cos(th)-mu*sp.sin(th))       # horizontal speed after landing (if >=0)
X = v**2*sp.sin(2*th)/g + u**2/(2*mu*g)
X2 = sp.simplify(sp.expand_trig(X))
print('X(theta) =', sp.simplify(sp.trigsimp(X)))
# stationary point
d = sp.simplify(sp.diff(X,th))
print('dX/dtheta at tan(theta)=mu :', sp.simplify(d.subs(th, sp.atan(mu))))
print('X at tan(theta)=mu:', sp.simplify(X.subs(th, sp.atan(mu))))
# alternative derivation: rewrite in 2theta
a = sp.symbols('a')
Xa = v**2/g*(sp.sin(a)/2 + sp.cos(a)*(1-mu**2)/(4*mu) + (1+mu**2)/(4*mu))
print('check identity:', sp.simplify(Xa.subs(a,2*th) - X))
# numeric brute force over theta, piecewise
def Xn(t,muv,vv=1.0,gv=1.0):
    uu = vv*(np.cos(t)-muv*np.sin(t))
    return vv**2*np.sin(2*t)/gv + (uu**2/(2*muv*gv) if uu>0 else 0.0)
ths=np.linspace(0,np.pi/2,200001)
for muv in [0.1,0.3,0.5,0.8,1.0,1.5,3.0]:
    uu=np.clip(np.cos(ths)-muv*np.sin(ths),0,None); vals=np.sin(2*ths)+uu**2/(2*muv); i=vals.argmax()
    pred = (1+muv**2)/(2*muv) if muv<1 else 1.0
    pth = np.degrees(np.arctan(muv)) if muv<1 else 45.0
    print(f'mu={muv}: numeric max {vals[i]:.6f} at {np.degrees(ths[i]):.3f} deg ; predicted {pred:.6f} at {pth:.3f} deg')
# impulse-level check: simulate collision with large normal force profile (e=0) and coulomb friction
def collide(theta,muv,vv=1.0):
    vx,vy=vv*np.cos(theta),-vv*np.sin(theta)  # landing velocity
    dt=1e-4; K=1e2  # N(t)=K*(-vy): damper contact, no rebound; gravity neglected during contact
    while vy<-1e-10:
        N=K*(-vy)
        vy+= N*dt
        if vx>0:
            vx=max(vx-muv*N*dt,0.0)
    return vx
for th_deg,muv in [(20,0.3),(40,0.5),(60,0.8),(30,2.0)]:
    t=np.radians(th_deg); print('collision sim',th_deg,muv,collide(t,muv),' formula',max(np.cos(t)-muv*np.sin(t),0))
