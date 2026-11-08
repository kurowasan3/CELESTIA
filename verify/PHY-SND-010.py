import numpy as np
from scipy.optimize import minimize_scalar, brentq
c=340.0; cR=350.0; L=10.0; D=5.0; f=500.0; W=1.7
def tau(y):
    r=minimize_scalar(lambda b: np.hypot(L,b)/c+np.hypot(D,y-b)/cR, bounds=(-1e3,1e3), method='bounded', options={'xatol':1e-12})
    return r.fun, r.x
for y in (2.0,8.775,30.0,300.0):
    t,b=tau(y); h=1e-4
    dtdy=(tau(y+h)[0]-tau(y-h)[0])/(2*h)
    fQ=f*(1+W*dtdy)   # Q moves -y: received phase rate
    sth=b/np.hypot(L,b); sph=(y-b)/np.hypot(D,y-b)
    print(f'y={y}: fQ={fQ:.5f}, f(1+W sin th/c)={f*(1+W*sth/c):.5f}, f(1+W sin ph/cR)={f*(1+W*sph/cR):.5f}, snell {sth/c:.6e} {sph/cR:.6e}')
print('limit fmax', f*(1+W/cR), ' critical sin', c/cR)
# position where sin th = 1/2
th=np.radians(30); b=L*np.tan(th); sph=cR*0.5/c; y=b+D*np.tan(np.arcsin(sph)); print('y at sin th=1/2:',y, ' fQ', f*(1+W/(2*c)))
# inverse: fmax=502.40 -> cR
print('cR from fmax 502.40:', W*f/(502.40-f))
