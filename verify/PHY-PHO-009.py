# T20213-1: ray-optics trap; small-angle force per ray, lens equivalence, 2-beam equilibrium under gravity
import sympy as sp, numpy as np
from scipy.optimize import fsolve
n,r,b,P,c=sp.symbols('n r b P c',positive=True)
th=sp.asin(b/r); ph=sp.asin(b/(n*r))
delta=2*(th-ph)
print('delta series:',sp.series(delta,b,0,3))       # 2(n-1)b/(n r)
f=n*r/(2*(n-1)); print('b/f =',sp.simplify(b/f))
# exact force per ray (vector): Delta p = p_out - p_in, magnitude 2 p sin(delta/2)
def ray_force(O, F0, u, Pw, nn, rr, cc):
    # beam line through F0 with unit dir u; sphere centre O. returns force on sphere
    w=O-F0; perp=w-np.dot(w,u)*u      # vector from line to O (perpendicular)
    bb=np.linalg.norm(perp)
    if bb<1e-30: return np.zeros(2)
    dl=2*(np.arcsin(bb/rr)-np.arcsin(bb/(nn*rr)))
    e=-perp/bb   # unit from O toward the line
    # outgoing direction rotated toward O (away from e side): p_out = cos dl u - sin dl e
    pin=Pw/cc*u; pout=Pw/cc*(np.cos(dl)*u-np.sin(dl)*e)
    return -(pout-pin)
nn=1.5; rr=2e-6; rho=1.05e3; Pw=10e-3; al=np.radians(30); beta=0.5; G=9.8; cc=3.0e8
M=4/3*np.pi*rr**3*rho; W=M*G
u1=np.array([np.sin(al),np.cos(al)]); u2=np.array([-np.sin(al),np.cos(al)])
F0=np.zeros(2)
def net(x):
    O=np.array(x)
    return ray_force(O,F0,u1,Pw,nn,rr,cc)+ray_force(O,F0,u2,beta*Pw,nn,rr,cc)+np.array([0,-W])
sol=fsolve(lambda z: net(np.array(z)*1e-6)/W,[0,-0.1],xtol=1e-13)*1e-6
k1=2*(nn-1)*Pw/(nn*cc*rr); k2=beta*k1
d1=W/(2*k1*np.sin(al)); d2=W/(2*k2*np.sin(al))
xa=(d1-d2)/(2*np.cos(al)); ya=-(d1+d2)/(2*np.sin(al))
print('M=',M,'Mg=',W,'k1=',k1)
print('numeric eq:',sol,' analytic:',xa,ya)
print('d1,d2 (um)=',d1*1e6,d2*1e6,' F per beam',W/(2*np.sin(al)))
# check each force magnitude equal
O=sol; print('|F1|,|F2|=',np.linalg.norm(ray_force(O,F0,u1,Pw,nn,rr,cc)),np.linalg.norm(ray_force(O,F0,u2,beta*Pw,nn,rr,cc)))
# stiffness ratio for equal beams: kx=2k cos^2, ky=2k sin^2 -> Tx/Ty=tan(alpha)
def Fx(dx,dy,bt=1.0):
    O=np.array([dx,dy]); return ray_force(O,F0,u1,Pw,nn,rr,cc)+ray_force(O,F0,u2,bt*Pw,nn,rr,cc)
h=1e-9
kx=-(Fx(h,0)[0]-Fx(-h,0)[0])/(2*h); ky=-(Fx(0,h)[1]-Fx(0,-h)[1])/(2*h)
print('kx/(2k cos^2)=',kx/(2*k1*np.cos(al)**2),' ky/(2k sin^2)=',ky/(2*k1*np.sin(al)**2),' Tx/Ty=',np.sqrt(ky/kx),' tan=',np.tan(al))
# --- with the parallel (2nd-order) component of Delta p removed, the analytic result should be reproduced
def ray_force_perp(O,F0,u,Pw,nn,rr,cc):
    Fv=ray_force(O,F0,u,Pw,nn,rr,cc); return Fv-np.dot(Fv,u)*u
def net2(x):
    O=np.array(x)
    return ray_force_perp(O,F0,u1,Pw,nn,rr,cc)+ray_force_perp(O,F0,u2,beta*Pw,nn,rr,cc)+np.array([0,-W])
sol2=fsolve(lambda z: net2(np.array(z)*1e-6)/W,[0,-0.1],xtol=1e-13)*1e-6
print('perp-only numeric eq:',sol2,' analytic:',xa,ya, ' rel.diff',(sol2-np.array([xa,ya]))/np.array([xa,ya]))
print('relative size of neglected axial part ~ delta/2 =',(nn-1)*d1/(nn*rr))
