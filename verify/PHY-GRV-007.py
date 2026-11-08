# T20221-1: contact force between two touching spheres in orbit: non-rotating (2) vs synchronously rotating (3)
import sympy as sp
G,M,m,s,r,phi,eps=sp.symbols('G M m s r phi epsilon',positive=True)
# planet at origin. pair CM at (r,0). outer sphere at CM + s*e, e=(cos phi, sin phi)
e=sp.Matrix([sp.cos(phi),sp.sin(phi)])
def grav(pos):
    d=sp.sqrt(pos.dot(pos)); return -G*M*pos/d**3
R0=sp.Matrix([r,0])
p1=R0+s*e; p2=R0-s*e
a_cm=(grav(p1)+grav(p2))/2
# case A: no rotation -> both spheres have acceleration a_cm. Outer sphere: m a_cm = m g(p1) - Gm^2/(4s^2) e + N e
NA=sp.simplify((m*a_cm-m*grav(p1)).dot(e)+G*m**2/(4*s**2))
NA1=sp.series(NA.subs(s,eps*s),eps,0,2).removeO().subs(eps,1)
print('N_A(phi) ~',sp.simplify(sp.expand(NA1)))
# case B: pair rotates with omega^2=GM/r^3 about planet keeping orientation phi relative to radial (co-rotating rigid)
w2=G*M/r**3
a1=-w2*p1   # uniform circular motion about origin for every point of the rigid co-rotating body
NB=sp.simplify((m*a1-m*grav(p1)).dot(e)+G*m**2/(4*s**2))
NB1=sp.series(NB.subs(s,eps*s),eps,0,2).removeO().subs(eps,1)
print('N_B(phi) ~',sp.simplify(sp.expand(NB1)))
print('N_B(0)=',sp.simplify(NB1.subs(phi,0)),'  N_B(pi/2)=',sp.simplify(NB1.subs(phi,sp.pi/2)))
# NB: exact CM orbit omega uses a_cm; first-order equal. Roche-type limit
rhoP,rho,Rp=sp.symbols('rho_P rho R_p',positive=True)
Ms=sp.Rational(4,3)*sp.pi*Rp**3*rhoP; ms=sp.Rational(4,3)*sp.pi*s**3*rho
rcB=sp.solve(sp.Eq(G*ms**2/(4*s**2),3*G*Ms*ms*s/r**3),r)
rcA=sp.solve(sp.Eq(G*ms**2/(4*s**2),2*G*Ms*ms*s/r**3),r)
print('rcB=',rcB,' rcA=',rcA)

# ---- (5) tidal torque on the pair about G and libration in synchronous rotation ----
eperp=sp.Matrix([-sp.sin(phi),sp.cos(phi)])
torque=(s*(m*grav(p1)).dot(eperp)) + ((-s)*(m*grav(p2)).dot(eperp))   # r x F (z) = arm * F.eperp
tq1=sp.series(torque.subs(s,eps*s),eps,0,3).removeO().subs(eps,1)
tq1=sp.simplify(sp.expand(tq1)); print('(5) torque ~',tq1)
assert sp.simplify(tq1-(-6*G*M*m*s**2/r**3*sp.sin(phi)*sp.cos(phi)))==0
I_=2*(sp.Rational(2,5)*m*s**2+m*s**2)
k=sp.simplify(-sp.diff(tq1,phi).subs(phi,0)/I_/(G*M/r**3)); print('(5) Omega_lib^2/omega^2 at phi=0:',k)
assert k==sp.Rational(15,7)
print('(5) stiffness at phi=pi/2:',sp.simplify(-sp.diff(tq1,phi).subs(phi,sp.pi/2)))
print('(5) period ratio T_lib/T =',sp.sqrt(1/k),float(sp.sqrt(1/k)))
# numeric: integrate full planar rigid-body + orbit (exact gravity), small libration, measure period
import numpy as np
from scipy.integrate import solve_ivp
Gv,Mv,mv,sv,rv=1.0,1.0,1e-9,1e-3,1.0
Iv=2*(0.4*mv*sv**2+mv*sv**2)
def acc(P): d=np.linalg.norm(P); return -Gv*Mv*P/d**3
def rhs(t,y):
    X,Y,VX,VY,psi,w=y; Rg=np.array([X,Y]); e=np.array([np.cos(psi),np.sin(psi)]); ep=np.array([-np.sin(psi),np.cos(psi)])
    F1=mv*acc(Rg+sv*e); F2=mv*acc(Rg-sv*e)
    a=(F1+F2)/(2*mv); tq=sv*F1.dot(ep)-sv*F2.dot(ep)
    return [VX,VY,a[0],a[1],w,tq/Iv]
om=np.sqrt(Gv*Mv/rv**3); phi0=0.02
y0=[rv,0,0,om*rv,phi0,om]
so=solve_ivp(rhs,[0,40],y0,rtol=1e-11,atol=1e-13,dense_output=True,max_step=0.01)
ts=np.linspace(0,40,40001); Y=so.sol(ts)
ph=Y[4]-np.unwrap(np.arctan2(Y[1],Y[0]))
zc=[ts[i] for i in range(1,len(ts)) if ph[i-1]>0>=ph[i]]
Tl=np.mean(np.diff(zc)); Torb=2*np.pi/om
print('(5) numeric T_lib/T_orb=%.4f theory=%.4f'%(Tl/Torb,np.sqrt(7/15)))
assert abs(Tl/Torb-np.sqrt(7/15))<2e-3
