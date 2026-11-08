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
print('numeric rho_P/rho=2:',float((24)**(sp.Rational(1,3))),float(16**(sp.Rational(1,3))))
