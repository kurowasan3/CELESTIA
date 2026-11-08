# K20043-1 verification: electromagnetic swing, vertical vs horizontal field damping
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
from scipy.signal import argrelextrema
# symbolic average <theta^2 theta'^2>
t,W,th0=sp.symbols('t Omega theta0',positive=True)
th=th0*sp.cos(W*t)
avg=sp.simplify(sp.integrate((th*sp.diff(th,t))**2,(t,0,2*sp.pi/W))/(2*sp.pi/W))
print('<theta^2 theta_dot^2> =',avg)
m,a,g,B,L,R=0.05,0.10,9.8,0.5,0.20,0.010   # gamma=(BL)^2/(mR)
K=(B*L*a)**2/R
W0=np.sqrt(g/a)
print('(BL)^2/(mR)=',(B*L)**2/(m*R),' Omega=',W0)
def run(kind,theta0,tend):
    def rhs(tt,y):
        q,w=y
        f=np.cos(q)**2 if kind=='v' else np.sin(q)**2
        return [w,(-m*g*a*np.sin(q)-K*f*w)/(m*a*a)]
    sol=solve_ivp(rhs,[0,tend],[theta0,0],rtol=1e-10,atol=1e-12,max_step=0.005)
    return sol
# vertical field: exponential decay of energy with tau=mR/(BL)^2
Bsave=B; B=0.05; K=(B*L*a)**2/R
th0v=0.05; sol=run('v',th0v,6.0)
E=0.5*m*a*a*sol.y[1]**2+m*g*a*(1-np.cos(sol.y[0]))
tau=m*R/(B*L)**2
i=np.searchsorted(sol.t,5.0)
print('vertical (B=0.05, weak damping): E(5s)/E0 num',E[i]/E[0],' exp(-t/tau)',np.exp(-sol.t[i]/tau))
B=Bsave; K=(B*L*a)**2/R
# horizontal: amplitude law theta0/sqrt(1+(BL theta0)^2 t/(4mR))
th0h=0.20; tend=200.0; sol=run('h',th0h,tend)
q=sol.y[0]; idx=argrelextrema(q,np.greater)[0]
for tt in (50,100,200):
    j=idx[np.argmin(abs(sol.t[idx]-tt))]
    pred=th0h/np.sqrt(1+(B*L*th0h)**2*sol.t[j]/(4*m*R))
    print(f'horizontal t={sol.t[j]:.1f}: amp num {q[j]:.5f} formula {pred:.5f}')
# half-amplitude times
print('t_half vertical =',2*np.log(2)*m*R/(B*L)**2,' horizontal(theta0=0.10)=',12*m*R/((B*L)**2*0.01),' ratio',6/(np.log(2)*0.01))
