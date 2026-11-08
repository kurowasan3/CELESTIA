# N23-1: RC with R2 parallel to C; determine R1,R2,C from graphs; half-times ratio
import numpy as np
from scipy.integrate import solve_ivp, quad
V=6.0; R1=2e3; R2=4e3; C=3e-3
f=lambda t,y:[((V-y[0]/C)/R1 - (y[0]/C)/R2)]
s=solve_ivp(f,[0,200],[0.0],max_step=0.01,rtol=1e-11,atol=1e-14,dense_output=True)
IC=lambda t:((V-s.sol(t)[0]/C)/R1-(s.sol(t)[0]/C)/R2)
print('IC(0) mA',IC(0)*1e3,'-> R1=V/IC0 =',V/IC(0))
qf=s.sol(200)[0]; Vf=qf/C; print('Vf',Vf,'qf mC',qf*1e3)
print('area to 20 s (squares of 0.2mC):',quad(IC,0,20,limit=500)[0]/0.2e-3, ' to inf:',quad(IC,0,200,limit=500)[0]/0.2e-3)
# discharge through R2
g=lambda t,y:[-(y[0]/C)/R2]
s2=solve_ivp(g,[0,300],[qf],max_step=0.01,rtol=1e-11,atol=1e-14,dense_output=True)
Id0=qf/C/R2; print('discharge I0 mA',Id0*1e3,' R1+R2=V/Id0=',V/Id0)
from scipy.optimize import brentq
tc=brentq(lambda t:IC(t)-IC(0)/2,0.01,50); td=brentq(lambda t:s2.sol(t)[0]/C/R2-Id0/2,0.01,100)
print('tc',tc,'td',td,'ratio',td/tc)
Vc=s.sol(tc)[0]/C; print('Vc at tc',Vc,'battery current mA',(IC(tc)+Vc/R2)*1e3)
