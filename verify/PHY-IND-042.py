# Q181-1 self-excited Faraday disk dynamo
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
B0,mu,n,r,R,w,I,tau,L,J=sp.symbols('B0 mu0 n r R omega I tau L J',real=True)
Isol=sp.solve(sp.Eq((B0+mu*n*I)*w*r**2/2,I*R),I)[0]; print('I(w)=',sp.simplify(Isol))
wc=2*R/(mu*n*r**2); print('wc',wc)
print('w->-oo I:',sp.limit(Isol.subs({B0:1,mu:1,n:1,r:1,R:1}),w,-sp.oo), ' field total:',sp.limit((B0+mu*n*Isol).subs({mu:1,n:1,r:1,R:1}),w,-sp.oo))
# constant torque, B0=0: tau = mu n I^2 r^2/2 ; emf: mu n I w r^2/2 = I R
sol=sp.solve([sp.Eq(tau,mu*n*I**2*r**2/2),sp.Eq(mu*n*I*w*r**2/2,I*R)],[I,w],dict=True); print(sol)
# energy check
Iv=sp.sqrt(2*tau/(mu*n*r**2)); print('power check',sp.simplify(tau*wc-Iv**2*R))
# with B0: torque tau = (B0+mu n I) I r^2/2, emf (B0+mu n I) w r^2/2 = I R  -> w as function of tau
# numeric dynamics: J dw/dt = tau - (B0+mu n I) I r^2/2 ; L dI/dt = (B0+mu n I) w r^2/2 - I R
p=dict(B0=5e-5,mu=4e-7*np.pi,n=1e4,r=0.1,R=0.05,L=0.02,J=0.05,tau=2.0)
mu_n=p['mu']*p['n']; wcv=2*p['R']/(mu_n*p['r']**2)
def f(t,y):
    om,Ii=y; Bt=p['B0']+mu_n*Ii
    return [(p['tau']-Bt*Ii*p['r']**2/2)/p['J'], (Bt*om*p['r']**2/2-Ii*p['R'])/p['L']]
so=solve_ivp(f,[0,400],[0.0,0.0],rtol=1e-10,atol=1e-12,method='LSODA')
print('final omega',so.y[0,-1],' wc',wcv,' final I',so.y[1,-1],' I(B0=0 theory)',np.sqrt(2*p['tau']/(mu_n*p['r']**2)))
# fixed omega below wc: I formula numeric
for om in [0.5*wcv,-0.5*wcv,-50*wcv]:
    print(om/wcv, p['B0']*om*p['r']**2/(2*p['R']-mu_n*om*p['r']**2), ' -B0/mu n=',-p['B0']/mu_n)
# --- dynamics: with constant torque and no friction the (omega,I) system oscillates around the steady state;
# time-average of omega over a long time equals wc (d lnI/dt = (mu n r^2/2L)(omega-wc) when B0->0)
p['B0']=1e-9
so=solve_ivp(f,[0,4000],[0.0,0.0],rtol=1e-11,atol=1e-14,method='LSODA',dense_output=True)
tt=np.linspace(1000,4000,600001); om=so.sol(tt)[0]; Ii=so.sol(tt)[1]
print('time-avg omega',np.mean(om),' wc',wcv,'  mean torque-current check <I^2>', np.mean(Ii**2), ' 2tau/(mu n r^2)=',2*p['tau']/(mu_n*p['r']**2))
# with viscous friction the motion converges to a steady state (here torque balance includes friction)
