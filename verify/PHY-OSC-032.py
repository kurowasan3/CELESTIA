# J034-1: Coulomb-friction spring oscillator, inverse determination
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,k,g,X0,xf,T=0.50,50.0,9.8,0.140,0.020,0.94
half=np.pi*np.sqrt(m/k); n=T/half; print('half period',half,' T/half',n)
n=3
d=(X0-xf)/(2*n); mup=k*d/(m*g); print('d',d,"mu'",mup)
mu_lo=k*xf/(m*g); mu_hi=k*(X0-4*d)/(m*g); print('mu range',mu_lo,mu_hi)
# path and heat
turns=[X0-2*j*d for j in range(n+1)]; L=sum(turns[j]+turns[j+1] for j in range(n)); print('turns',turns,'L',L)
print('heat',mup*m*g*L,' energy diff',0.5*k*(X0**2-xf**2))
# simulate with mu chosen in range, e.g. 0.4
def simulate(mu,mup,X0):
    x,v,t=X0,0.0,0.0; hist=[]
    while True:
        if abs(k*x)<=mu*m*g and v==0.0: return x,t,hist
        s=-np.sign(k*x) # direction of motion starting from rest
        def f(tt,y): return [y[1], (-k*y[0]-mup*m*g*s)/m]
        ev=lambda tt,y:y[1]; ev.terminal=True; ev.direction=-s
        sol=solve_ivp(f,[t,t+10],[x,1e-12*s],events=ev,rtol=1e-11,atol=1e-13,max_step=1e-3)
        x=sol.y_events[0][0][0]; t=sol.t_events[0][0]; v=0.0; hist.append((round(t,4),round(x,5)))
for mu in [0.21,0.4,0.6]:
    print('mu',mu,simulate(mu,mup,X0))
# (5) without timing: possible mu' = 3kX0/(7 n m g) for odd n, check n=1,5 by simulation with mu = mu' (min allowed)
for nn in [1,3,5,7]:
    dd=3*X0/(7*nn); mm=k*dd/(m*g); xf_,t_,h=simulate(max(mm,k*xf/(m*g))*1.0000001,mm,X0)
    print(nn,'mu_p',mm,'final',xf_,'t/half',t_/half)
