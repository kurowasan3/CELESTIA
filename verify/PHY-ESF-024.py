# S0503-1 verification: energies via (1/2)sum QV and independently via field energy integral eps0 E^2/2
import sympy as sp
Q,r,R,Rp,k,rho=sp.symbols('Q r R Rp k rho',positive=True)
eps0=1/(4*sp.pi*k)
def field_energy(pieces):  # pieces: list of (E(rho), lo, hi)
    return sp.simplify(sum(sp.integrate(eps0*E**2/2*4*sp.pi*rho**2,(rho,lo,hi)) for E,lo,hi in pieces))
Ec=k*Q/rho**2
U0=field_energy([(Ec,r,sp.oo)])                       # C alone
U1=field_energy([(Ec,r,R),(Ec,Rp,sp.oo)])             # neutral shell closed
U2=field_energy([(Ec,r,R)])                           # after grounding (outside field 0)
U3=field_energy([(Ec,r,sp.oo),(-k*Q/rho**2,Rp,sp.oo)])# C far away + shell (-Q) alone (separate, no cross term)
# (1/2) sum QV
V1C=k*Q*(1/r-1/R+1/Rp); U1b=sp.Rational(1,2)*Q*V1C
V2C=k*Q*(1/r-1/R); U2b=sp.Rational(1,2)*Q*V2C
U3b=k*Q**2/(2*r)+k*Q**2/(2*Rp)
for a_,b_ in [(U0,k*Q**2/(2*r)),(U1,U1b),(U2,U2b),(U3,U3b)]:
    assert sp.simplify(a_-b_)==0
W_elec=sp.simplify(U0-U1); J=sp.simplify(U1-U2); W3=sp.simplify(U3-U2)
print('work by electric force when closing:',sp.factor(W_elec))
print('Joule heat on grounding:',J)
print('external work to remove C:',sp.factor(W3))
assert sp.simplify(W_elec-k*Q**2/2*(1/R-1/Rp))==0
assert sp.simplify(J-k*Q**2/(2*Rp))==0
assert sp.simplify(W3-k*Q**2/2*(1/R+1/Rp))==0
# numeric sanity (field energy integral numerically)
import numpy as np
from scipy.integrate import quad
kk,QQ,rr,RR,RRp=9e9,1e-8,0.01,0.05,0.08
f=lambda x: (1/(8*np.pi*kk))*(kk*QQ/x**2)**2*4*np.pi*x**2
num=quad(f,RR,RRp)[0]; print('numeric closing work',num,' formula',kk*QQ**2/2*(1/RR-1/RRp))
print('ALL OK')
