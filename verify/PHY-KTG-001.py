# M48-1: photon gas in a mirrored sphere -> P = E/(3V); compare monatomic gas
import sympy as sp, numpy as np
eps,c,r,th,t,E,N=sp.symbols('epsilon c r theta t E N',positive=True)
imp=2*(eps/c)*sp.cos(th); n=c*t/(2*r*sp.cos(th))
F1=sp.simplify(imp*n/t); print('force per photon',F1)
P=sp.simplify(E/r/(4*sp.pi*r**2)); V=sp.Rational(4,3)*sp.pi*r**3
print('P*V/E',sp.simplify(P*V/E))
# general particle momentum p speed v: force per particle p v / r ; massive p v = m v^2 = 2K
# ray tracing check
rng=np.random.default_rng(0)
R=1.0; cc=1.0; Ttot=200.0
def impulse_rate(pos,d):
    tt=0; J=0
    while True:
        b=pos@d; q=pos@pos-R*R
        s=-b+np.sqrt(max(b*b-q,0.0))
        if tt+s>Ttot: break
        tt+=s; pos=pos+s*d
        nrm=pos/np.linalg.norm(pos); vn=d@nrm
        J+=2*abs(vn)  # eps/c=1
        d=d-2*vn*nrm
    return J/Ttot
vals=[]
for i in range(200):
    p=rng.normal(size=3); p=p/np.linalg.norm(p)*R*rng.random()**(1/3)
    d=rng.normal(size=3); d/=np.linalg.norm(d)
    vals.append(impulse_rate(p,d))
print('mean radial impulse rate per photon (expect eps/r=1):',np.mean(vals),np.std(vals))
# numeric (5): E=1.0 J, r=5.0 cm
rn=0.050; En=1.0
print('P =',En/(3*4/3*np.pi*rn**3),'Pa')
