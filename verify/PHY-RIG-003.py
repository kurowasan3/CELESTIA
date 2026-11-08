# M06-1: disk with off-center hole on rough incline
import sympy as sp, numpy as np, math
R,M0,M,g,al,ph,m=sp.symbols('R M0 M g alpha phi m',positive=True)
# (1) CG by negative mass: hole radius R/2 at (R/3,0)
mh=M0*sp.Rational(1,4); xG=sp.simplify((0*M0 - mh*R/3)/(M0-mh)); print('xG',xG)
# numeric Monte Carlo CG check
rng=np.random.default_rng(0); pts=rng.uniform(-1,1,(4_000_000,2)); r=1
inside=(pts[:,0]**2+pts[:,1]**2<1)&((pts[:,0]-1/3)**2+pts[:,1]**2>0.25)
print('MC xG',pts[inside,0].mean(),' pred',-1/9)
e=R/9
# (2) potential: U/Mg = -R*phi*sin(al) - e*cos(psi0+phi); equilibrium & stability
psi=sp.symbols('psi')
U=-R*ph*sp.sin(al)-e*sp.cos(psi+ph)
dU=sp.diff(U,ph).subs(ph,0); d2U=sp.diff(U,ph,2).subs(ph,0)
print('dU=0:',sp.solve(sp.Eq(dU,0),sp.sin(psi)) if False else sp.simplify(dU), ' d2U',d2U)
# brute-force numeric: rigid body torque about contact point for sin(al)=1/18 -> sin psi=1/2
a=math.asin(1/18); ev=1/9
for psi_v in [math.pi/6, 5*math.pi/6]:
    # O at origin, contact P = -(sin a, cos a) (slope descending to right). G = e*( -sin psi, -cos psi)  (psi from downward vertical toward uphill/left)
    Gx=-ev*math.sin(psi_v); Px=-math.sin(a)
    print('psi',psi_v,' G_x-P_x =',Gx-Px, ' second deriv sign', ev*math.cos(psi_v))
# friction: contact force must equal Mg vertical -> f = Mg sin a, N = Mg cos a
print('mu >= tan alpha (contact force vertical)')
# (3) mass m at rim on ray OG, need offset >= R/3
off=(M*R/9+m*R)/(M+m); print('m min',sp.solve(sp.Eq(off,R/3),m))
# check optimality by random placement
best=0
for k in range(200000):
    th=rng.uniform(0,2*math.pi); rr=math.sqrt(rng.uniform(0,1))
    x,y=rr*math.cos(th),rr*math.sin(th)
    if (x-1/3)**2+y**2<0.25: continue
    mm=1/3; gx=(-1/9+mm*x)/(1+mm); gy=mm*y/(1+mm); best=max(best,math.hypot(gx,gy))
print('max offset with m=M/3 random placements',best,' (<=1/3)')
