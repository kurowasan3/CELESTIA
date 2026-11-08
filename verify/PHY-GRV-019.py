# Q111-1: solid hemisphere: field at center of flat face and at the pole
import numpy as np, sympy as sp
from scipy.integrate import tplquad
G=rho=R=1.0
g0=4*np.pi/3
# field z-component at point (0,0,z0) from hemisphere z>=0, r<=R (spherical coords)
def gz(z0):
    f=lambda r,th,ph: r**2*np.sin(th)*(r*np.cos(th)-z0)/((r*np.sin(th))**2+(r*np.cos(th)-z0)**2)**1.5
    val,err=tplquad(f,0,2*np.pi,0,np.pi/2,0,R,epsabs=1e-9,epsrel=1e-9)
    return val
gO=gz(0.0); print('g at O',gO,' theory pi =',np.pi,' ratio',gO/g0)
gP=-gz(1.0-1e-9); print('g at pole',gP,' theory',2/3*(3-np.sqrt(2))*np.pi,' ratio',gP/g0,(3-np.sqrt(2))/2)
# sympy: shell result and disk-sum result
th,s=sp.symbols('theta s',positive=True); Gs,sig,Rs=sp.symbols('G sigma R',positive=True)
shell=sp.integrate(Gs*sig*2*sp.pi*Rs**2*sp.sin(th)*sp.cos(th)/Rs**2,(th,0,sp.pi/2)); print('hemishell',shell)
disk=sp.integrate(2*sp.pi*Gs*(1-sp.sqrt(s/(2*Rs))),(s,0,Rs)); print('pole via disks',sp.simplify(disk))
z=sp.symbols('z',positive=True)
other=sp.integrate(2*sp.pi*Gs*(1-(Rs+z)/sp.sqrt(2*Rs*(Rs+z))),(z,0,Rs)); print('other hemisphere at P',sp.simplify(other))
print('superposition',sp.simplify(sp.Rational(4,3)*sp.pi*Gs*Rs-other - disk))
