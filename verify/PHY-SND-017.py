# O20173-1 tuning-fork model: 4 point sources (quadrupole) far-field amplitude
import numpy as np, sympy as sp
def amp(ka,phi,block_S=False):
    n=np.array([np.cos(phi),np.sin(phi)])
    src=[((0,1),1),((0,-1),1),((1,0),-1)]+([] if block_S else [((-1,0),-1)])
    return abs(sum(s*np.exp(1j*ka*np.dot(p,n)) for p,s in src))
# (2) formula check
for ka in [0.3,1.0,2.5]:
    for phi in np.linspace(0,2*np.pi,7):
        f=2*abs(np.cos(ka*np.sin(phi))-np.cos(ka*np.cos(phi)))
        assert abs(f-amp(ka,phi))<1e-12
print("formula 2A|cos(ka sin)-cos(ka cos)| OK")
# count zeros in [0,2pi) for a/lambda below/above 1/sqrt2
def zeros(a_over_lam):
    ka=2*np.pi*a_over_lam
    ph=np.linspace(0,2*np.pi,400001)[:-1]
    g=np.cos(ka*np.sin(ph))-np.cos(ka*np.cos(ph))
    return int(np.sum(np.sign(g)!=np.sign(np.roll(g,-1))))
for al in [0.3,0.6,0.70,0.72,0.9]:
    print("a/lambda=",al,"zero directions:",zeros(al))
# (3) small ka
ka=0.05
for phi in [0,0.3,np.pi/4,1.2]:
    print(phi, amp(ka,phi), ka**2*abs(np.cos(2*phi)))
x,p=sp.symbols('x phi')
print(sp.series(2*(sp.cos(x*sp.sin(p))-sp.cos(x*sp.cos(p))),x,0,4))
# (5) blocking S, phi=0
lam=340/425; a=0.02; ka=2*np.pi*a/lam
r=(amp(ka,0,True)/amp(ka,0))**2
print("ka=",ka," intensity ratio exact",r," approx (ka)^-4",ka**-4)
