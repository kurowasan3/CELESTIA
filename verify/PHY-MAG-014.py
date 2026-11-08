# T20231-1: decay of a charged particle in uniform B: midpoint of daughters' gyration centres = parent's centre
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
q,B,m=1.0,1.0,1.0; R0=1.0
V0=q*B*R0/(2*m)   # X: mass 4m, charge 2q: R0=4m V0/(2qB)
def centre(r,p,e):  # positive charge, B along +z: c = r + p x z /(eB)
    return r+np.array([p[1],-p[0]])/(e*B)
def traj(r,v,mass,e,T):
    def f(t,y): return [y[2],y[3], e/mass*(y[3]*B), e/mass*(-y[2]*B)]
    return solve_ivp(f,[0,T],[*r,*v],rtol=1e-11,atol=1e-12,dense_output=True)
# parent check: start at P=(0,R0) moving +x -> clockwise about O
P=np.array([0.0,R0]); t_hat=np.array([1.0,0.0])
print('parent centre',centre(P,4*m*V0*t_hat,2*q))
rng=np.random.default_rng(1)
for trial in range(3):
    ang=rng.uniform(0,2*np.pi); u=rng.uniform(0.2,3)*V0
    e_=np.array([np.cos(ang),np.sin(ang)])
    vA=V0*t_hat+u*e_; vB=V0*t_hat-u/3*e_
    cA=centre(P,m*vA,q); cB=centre(P,3*m*vB,q)
    print('random emission: (cA+cB)/2 =',(cA+cB)/2)
# forward emission: radii & tangency & meeting time
for u in [0.5*V0,2*V0,4*V0]:
    vA=(V0+u)*t_hat; vB=(V0-u/3)*t_hat
    w=m*u/(q*B)
    cA=centre(P,m*vA,q); cB=centre(P,3*m*vB,q)
    print('u=',u,' rA',np.linalg.norm(P-cA),' formula',abs(R0/2+w),' rB',np.linalg.norm(P-cB),' formula',abs(1.5*R0-w))
    TA=2*np.pi*m/(q*B)
    sA=traj(P,vA,m,q,3.2*TA); sB=traj(P,vB,3*m,q,3.2*TA)
    ts=np.linspace(0.01,3.2*TA,200001)
    dist=np.linalg.norm(sA.sol(ts)[:2]-sB.sol(ts)[:2],axis=0)
    i=np.argmin(dist); print('   closest approach t/TA=',ts[i]/TA,' dist',dist[i],' pos',sA.sol(ts[i])[:2])
    # check no earlier near-meeting
    early=dist[(ts>0.1*TA)&(ts<2.9*TA)].min(); print('   min dist for 0.1TA<t<2.9TA:',early)
# B passes through O when w=R0 -> u=2V0, Q=(2/3) m u^2
u=2*V0; vB=(V0-u/3)*t_hat; cB=centre(P,3*m*vB,q); rB=np.linalg.norm(P-cB)
print('u=2V0: |cB|=',np.linalg.norm(cB),' rB=',rB, ' Q=',2/3*m*u**2,' (2/3)q^2B^2R0^2/m=',2/3*q**2*B**2*R0**2/m)
# energy release check: in X frame, A speed u, B speed u/3: KE = (1/2)m u^2+(1/2)3m(u/3)^2
uu,mm=sp.symbols('u m',positive=True); print('Q(u)=',sp.simplify(mm*uu**2/2+3*mm*(uu/3)**2/2))
# special case u=V0: both daughter circles coincide with X's circle -> A laps B at t=1.5 TA at the opposite point
u=V0; vA=(V0+u)*t_hat; vB=(V0-u/3)*t_hat; TA=2*np.pi*m/(q*B)
sA=traj(P,vA,m,q,3.2*TA); sB=traj(P,vB,3*m,q,3.2*TA)
ts=np.linspace(0.1*TA,3.2*TA,200001); dist=np.linalg.norm(sA.sol(ts)[:2]-sB.sol(ts)[:2],axis=0)
i=np.argmin(dist); print('u=V0: first meeting t/TA=',ts[np.argmax(dist<1e-3)]/TA,' at',sA.sol(ts[np.argmax(dist<1e-3)])[:2])
