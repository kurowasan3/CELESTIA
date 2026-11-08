import numpy as np
from scipy.integrate import solve_ivp
# D1: P1 (0,0) lower, Q1 (0,2a) upper, vel +x ; D2: Q2 (0,0) upper, P2 (0,-4a) lower, vel -x
a=1.0; v=1.0; m=1.0
kr=4e5; kc=4e6; rad=0.01
L1,L2=2*a,4*a
def rhs(t,s):
    p=s[:8].reshape(4,2); u=s[8:].reshape(4,2); F=np.zeros((4,2))
    for (i,j,L) in [(0,1,L1),(2,3,L2)]:
        dvec=p[j]-p[i]; dist=np.linalg.norm(dvec); f=kr*(dist-L)*dvec/dist
        F[i]+=f; F[j]-=f
    for i in (0,1):
        for j in (2,3):
            dvec=p[j]-p[i]; dist=np.linalg.norm(dvec)
            if dist<2*rad:
                f=kc*(2*rad-dist)*dvec/dist; F[j]+=f; F[i]-=f
    return np.concatenate([u.ravel(),(F/m).ravel()])
# start slightly before contact so that P1 and Q2 approach head on
x0=0.05
p0=np.array([[-x0-rad,0],[-x0-rad,2*a],[x0+rad,0],[x0+rad,-4*a]])
u0=np.array([[v,0],[v,0],[-v,0],[-v,0]])
s0=np.concatenate([p0.ravel(),u0.ravel()])
T=10.0
sol=solve_ivp(rhs,[0,T],s0,method='DOP853',rtol=1e-9,atol=1e-11,max_step=2e-3,dense_output=True)
ts=np.linspace(0,T,2001)
for tt in [0.02,0.2,1.0,3.0,6.0,6.4,7.0,9.0]:
    s=sol.sol(tt); u=s[8:].reshape(4,2)
    c1=(u[0]+u[1])/2; c2=(u[2]+u[3])/2
    print(f't={tt:5.2f} vcm1={np.linalg.norm(c1):.4f} vcm2={np.linalg.norm(c2):.4f}', np.round(u,3).tolist())
print('expected 2nd collision at t ~',0.05/1+2*np.pi*a/v)
