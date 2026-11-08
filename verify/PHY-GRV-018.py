# Q110-1: spherical cavity -> uniform field; chord (Galileo) isochronism
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
rng=np.random.default_rng(0)
G=1.0; rho=1.0; R=1.0; d=0.45; b=0.4
w2=4*np.pi*G*rho/3
# (1) brute-force field at points in cavity by Monte Carlo over the solid (sphere minus cavity)
N=4_000_000
pts=rng.uniform(-R,R,(N,3)); pts=pts[(np.linalg.norm(pts,axis=1)<R)]
Op=np.array([d,0,0])            # cavity center on +x
pts=pts[np.linalg.norm(pts-Op,axis=1)>b]
dm=rho*(8*R**3)/N
for P in [Op, Op+np.array([0.2,0.1,-0.15]), Op+np.array([-0.3,0.,0.2])]:
    r=pts-P; E=G*dm*np.sum(r/np.linalg.norm(r,axis=1)[:,None]**3,axis=0)
    print('MC field',E,' theory',-w2*Op)
# (2) times along chords from P (cavity point farthest from O) to Q
gp=w2*d
P=Op+np.array([b,0,0])
def field(X): return -w2*X + w2*(X-Op)   # inside cavity
for phi in np.radians([0,20,45,70,85]):
    Q=Op+b*np.array([np.cos(np.pi-2*phi), np.sin(np.pi-2*phi),0])  # inscribed angle phi from diameter
    e=(Q-P)/np.linalg.norm(Q-P); Lc=np.linalg.norm(Q-P)
    def rhs(t,y): return [y[1], np.dot(field(P+y[0]*e),e)]
    ev=lambda t,y:y[0]-Lc; ev.terminal=True
    s=solve_ivp(rhs,[0,50],[0,0],events=ev,rtol=1e-11,atol=1e-12)
    print(f'phi={np.degrees(phi):.0f} angleQPO={np.degrees(np.arccos(-e[0])):.1f} T={s.t_events[0][0]:.6f} theory {2*np.sqrt(b/gp):.6f} vQ={s.y_events[0][0][1]:.5f} theory {gp*2*np.sqrt(b/gp)*np.cos(phi):.5f}')
# (4) b/R with b+d=R making T equal to quarter period of uniform sphere
bb=sp.symbols('b',positive=True); Rs=sp.symbols('R',positive=True)
sol=sp.solve(sp.Eq(4*bb/(Rs-bb),sp.pi**2/4),bb); print('b/R=',sp.simplify(sol[0]/Rs),float(sol[0].subs(Rs,1)))
