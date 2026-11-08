import numpy as np
A=2.0; T=0.40; w=2*np.pi/T
def check(direction,lam):
    k=2*np.pi/lam
    y=lambda x,t: A*np.sin(w*t-direction*k*x)
    return y,k
# graphs: P: A sin(wt), Q (x=3): -A cos(wt)
cands=[]
for d in (+1,-1):
    for lam in np.arange(3.01,100,0.001):
        k=2*np.pi/lam
        if abs(np.sin(-d*k*3)-(-1))<1e-7: cands.append((d,round(lam,3)))
print('candidates (dir,lambda) with lambda>3:',sorted(set(cands)))
# density excess  ~ -dy/dx at P, t=0
for d,lam in [(1,12.0),(-1,4.0)]:
    k=2*np.pi/lam; dydx=-d*k*A*np.cos(0); print('dir',d,'lam',lam,'density excess at P', -dydx>0 and '密' or '疎')
lam=12.0;k=2*np.pi/lam;v=lam/T; print('v',v,'cm/s')
xs=np.linspace(0,12,12001)
dens=lambda x,t: -(-k*A*np.cos(w*t-k*x))  # -dy/dx
print('t=0 most sparse x:',xs[np.argmin(dens(xs,0))])
vel=lambda x,t: A*w*np.cos(w*t-k*x)
i=np.argmin(vel(xs[:-1],0.10)); print('t=0.10 most -x velocity x:',xs[i],'density there',dens(xs[i],0.10))
ts=np.linspace(0,0.4,40001); print('Q first densest t:',ts[np.argmax(dens(3.0,ts))])
# relation dy/dt = -v dy/dx
print('max |dy/dt + v dy/dx|',np.max(np.abs(vel(xs,0.13)-v*dens(xs,0.13))))
assert abs(xs[i]-9)<1e-3 and abs(ts[np.argmax(dens(3.0,ts))]-0.10)<1e-4; print('OK')
