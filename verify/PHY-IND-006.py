# N41-1: ring of two semicircular resistors around a solenoid-like flux region; path-dependent voltmeter
import sympy as sp, numpy as np
k,b,R1,R2,r,rv=sp.symbols('k b R1 R2 r r_V',positive=True)
eps=k*sp.pi*b**2
I=eps/(R1+R2)  # clockwise ring current
# meshes, clockwise currents J1 (upper half-disk), J2 (lower half-disk); diameter element resistance x
def mesh(x):
    J1,J2=sp.symbols('J1 J2')
    sol=sp.solve([sp.Eq(J1*R1+(J1-J2)*x,eps/2),sp.Eq(J2*R2+(J2-J1)*x,eps/2)],[J1,J2])
    return sol[J1],sol[J2]
J1,J2=mesh(r)
iAC=sp.simplify(J2-J1)  # current along diameter from A to C
print('i_AC (resistor r) =',sp.factor(iAC))
J1v,J2v=mesh(rv); Vb=sp.limit((J2v-J1v)*rv,rv,sp.oo)
print('voltmeter on diameter reading (A->C current * rV) =',sp.factor(Vb), ' expected',sp.factor(eps*(R1-R2)/(2*(R1+R2))))
# outside-routes: loops with no flux
print('upper route  V =',sp.simplify(I*R1),'  lower route V =',sp.simplify(-I*R2))
print('difference upper-lower =',sp.simplify(I*R1+I*R2-eps))
# r -> 0 check
print('r->0 J1,J2:',[sp.simplify(sp.limit(e,r,0)) for e in mesh(r)])
# Independent numeric check: discretized circuit via line integral of induced E (E_phi = -k b^2/(2 rho) outside b)
# EMF along a path = -∫E·dl ; diameter path: E_phi perpendicular? on the diameter through center E is tangential -> perpendicular to diameter, so ∫E·dl along diameter = 0.
# Upper arc from A(-a,0) clockwise to C(a,0): ∫E·dl where E clockwise-direction magnitude k b^2/(2a): path length pi a
a=2.0; kv=1.0; bv=1.0; R1v=3.0; R2v=1.0
Ecw=kv*bv**2/(2*a)   # induced E along clockwise direction (flux out of page increasing)
emf_upper=Ecw*np.pi*a; emf_lower=Ecw*np.pi*a
# with E-field picture: in arc, IR1 = emf_upper + (phi_A - phi_C) ; ring current I=(emf_u+emf_l)/(R1+R2)
Iv=(emf_upper+emf_lower)/(R1v+R2v)
phiA_minus_phiC=Iv*R1v-emf_upper   # electrostatic potential difference
# voltmeter on diameter: no induced E along diameter, so reads electrostatic difference
print('E-field method diameter reading',phiA_minus_phiC,' formula',float((eps*(R1-R2)/(2*(R1+R2))).subs({k:kv,b:bv,R1:R1v,R2:R2v})))
