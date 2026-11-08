# N65-1: Coulomb barrier, CM energy, thermal mean relative KE, 2->1 impossibility
import sympy as sp, numpy as np
k0,e,R,kB=9.0e9,1.6e-19,4.0e-15,1.4e-23
U=k0*e**2/R; MeV=1.6e-13; print('U=',U,'J =',U/MeV,'MeV')
# head-on: projectile m1 on target m2 at rest; closest approach when common velocity
m1,m2,v,Uc=sp.symbols('m1 m2 v U',positive=True)
V=m1*v/(m1+m2); lost=sp.simplify(m1*v**2/2-(m1+m2)*V**2/2)
Kmin=sp.solve(sp.Eq(lost,Uc),v)[0]; Kp=sp.simplify(m1*Kmin**2/2); print('K_min=',Kp)
m=1
print('d on t:',U/MeV*5/3,' t on d:',U/MeV*5/2)
# numeric check by integrating Coulomb repulsion of point charges (d on t), K slightly above/below
from scipy.integrate import solve_ivp
mu=1.66e-27
def run(K,ma,mb):
    va=np.sqrt(2*K/ma)
    def f(t,y):
        xa,xb,ua,ub=y; F=k0*e**2/(xb-xa)**2
        return [ua,ub,-F/ma,F/mb]
    ev=lambda t,y:(y[1]-y[0])-R; ev.terminal=True
    turn=lambda t,y:(y[3]-y[2]); turn.terminal=True
    s=solve_ivp(f,[0,1e-18],[0,2e-12,va,0],events=[ev,turn],rtol=1e-10,atol=1e-22,max_step=1e-22)
    return 'touch' if len(s.t_events[0]) else 'turn'
for K in [U*5/3*0.99,U*5/3*1.01]:
    print('d on t, K/U=',K/U,run(K,2*mu,3*mu))
for K in [U*5/2*0.99,U*5/2*1.01]:
    print('t on d, K/U=',K/U,run(K,3*mu,2*mu))
# thermal: mean of mu*vrel^2/2 = 3/2 kT (Monte Carlo)
rng=np.random.default_rng(0); T=1.0; ma,mb=2.0,3.0
v1=rng.normal(0,np.sqrt(T/ma),(400000,3)); v2=rng.normal(0,np.sqrt(T/mb),(400000,3))
muu=ma*mb/(ma+mb); print('MC mean rel KE /(1.5kT)=',np.mean(0.5*muu*np.sum((v1-v2)**2,1))/1.5)
Tn=2*U/(3*kB); print('T needed',Tn,' naive',U/(3*kB))
