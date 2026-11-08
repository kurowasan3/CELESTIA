# N43-1: solenoid end field = half (superposition), mutual inductance via reciprocity, charge independent of L
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, quad
mu0=4e-7*np.pi; n=1000.; I=1.; a=0.02; l=1.0   # length 2l, l>>a
def Baxis(z):  # z measured from centre, finite solenoid on |z|<=l
    f=lambda zz:(zz)/np.sqrt(zz**2+a**2)
    return mu0*n*I/2*(f(z+l)-f(z-l))
print('B(centre)/(mu0 n I)=',Baxis(0)/(mu0*n*I),'  B(end)/(mu0 n I)=',Baxis(l)/(mu0*n*I))
# small coil Y: N turns, area s at centre/end: M = N s B/I
N=50; s=1e-5
MP=N*s*Baxis(0)/I; ME=N*s*Baxis(l)/I
print('M_P/(mu0 n N s)=',MP/(mu0*n*N*s),' M_E/(mu0 n N s)=',ME/(mu0*n*N*s))
# reciprocity check: flux linkage of solenoid X due to a small dipole (Y) at centre, computed directly:
# dipole moment m = N s iY; flux through a loop of radius a at axial distance z from dipole:
# Phi(z) = mu0 m a^2 /(2 (a^2+z^2)^{3/2})  (vector potential result)
iY=1.0; mdip=N*s*iY
Phi=lambda z: mu0*mdip*a**2/(2*(a**2+z**2)**1.5)
for zc,name,Mx in [(0.0,'P',MP),(l,'E',ME)]:
    link=quad(lambda zp: n*Phi(zp-zc),-l,l,limit=400,points=[zc])[0]
    print(name,'direct linkage/iY =',link/iY,'  M from B(axis) =',Mx)
# charge through R with self inductance L when Y (current i0) pulled from P to far away
L=0.05; R=2.0; i0=3.0
def Mt(t):  # some smooth pulling profile from P (t=0) outward over 0.01 s
    zc=l*3*(t/0.01)**2 if t<0.01 else 3*l
    return N*s*Baxis(zc)/1.0 if True else 0
def rhs(t,y):
    h=1e-7; dM=(Mt(t+h)-Mt(t-h))/(2*h)
    I_,q=y; return [(-dM*i0-R*I_)/L, I_]
sol=solve_ivp(rhs,[0,0.5],[0,0],max_step=1e-5,rtol=1e-9,atol=1e-14)
print('charge num',sol.y[1,-1],'  expected (M_P - M_far) i0/R =',(MP-Mt(1.0))*i0/R,' ~ mu0 n N s i0/R =',mu0*n*N*s*i0/R)
