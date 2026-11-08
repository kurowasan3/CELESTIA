# J084-1: moving bat, echoes from two stationary surfaces separated by d -> cancellation condition; wind case
import numpy as np
from scipy.optimize import brentq
V=340.; d=0.004; v=8.0
def te_of(t,s,vb,x0=-5.0):
    xb=lambda tt: x0+vb*tt
    g=lambda te: V*(t-te)-((s-xb(te))+(s-xb(t)))
    return brentq(g,-10,t)
for vb in (v,-v):
    t=0.2
    dt=te_of(t,0.0,vb)-te_of(t,d,vb)
    # destructive if f*dt = 1/2 -> f = 1/(2 dt)
    f=1/(2*dt); print('vb',vb,' f_min_cancel',f,' formula',(V-vb)/(4*d))
    # echo frequency heard by bat
    h=1e-6; fe=f*(te_of(t+h,0.0,vb)-te_of(t-h,0.0,vb))/(2*h); print('   heard',fe,' formula',(V+vb)/(4*d))
f1=(V-v)/(4*d); f2=(V+v)/(4*d)
print('inverse d',V/(2*(f1+f2)),' v',V*(f2-f1)/(f1+f2))
# wind w toward cliff, bat at rest: phase difference f*(d/(V+w)+d/(V-w))
w=10.
print('wind f_cancel',1/(2*(d/(V+w)+d/(V-w))),' formula',(V**2-w**2)/(4*d*V))
print('numbers: d=V/(2(f1+f2)) with f1=20.75k? ', 340/(2*(41.5e3+43.5e3)), 340*(43.5-41.5)/(85))
