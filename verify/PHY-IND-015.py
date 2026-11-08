# T20012-1 verification: conducting ring falling along the axis of Bz = B0 + C z
import sympy as sp
import numpy as np

r, z, B0, C, a, m, g, R, v = sp.symbols('r z B0 C a m g R v', real=True)
Br = sp.Function('Br')
# div B = (1/r) d(r Br)/dr + dBz/dz = 0 with Br(0)=0
Bz = B0 + C*z
Br_sol = sp.integrate(-r*sp.diff(Bz, z), (r, 0, r))/r
print('B_r =', sp.simplify(Br_sol))
# flux & emf: ring moving down with speed v (z' = -v)
Phi = sp.pi*a**2*Bz
emf = sp.simplify(-sp.diff(Phi, z)*(-v))   # -dPhi/dt, positive = ccw from above
I = emf/R
print('emf (ccw positive) =', emf, '  I =', I)
# force on ring by direct integration of I dl x B, dl = a dphi phihat (ccw)
ph = sp.symbols('phi')
rhat = sp.Matrix([sp.cos(ph), sp.sin(ph), 0]); phat = sp.Matrix([-sp.sin(ph), sp.cos(ph), 0]); zhat = sp.Matrix([0, 0, 1])
Bvec = Br_sol.subs(r, a)*rhat + Bz*zhat
dF = I*a*phat.cross(Bvec)
F = sp.simplify(sp.Matrix([sp.integrate(dF[i], (ph, 0, 2*sp.pi)) for i in range(3)]))
print('net force =', F.T)
print('power check F_z v == I^2 R :', sp.simplify(F[2]*v - I**2*R) == 0)
vt = sp.solve(sp.Eq(F[2], m*g), v)
print('terminal speed =', vt)
# hoop: radial force per length = I*(phat x Bz zhat) . rhat = I Bz ; tension T = f * a
f_rad = sp.simplify((phat.cross(Bz*zhat)).dot(rhat))*I
T = sp.simplify((f_rad*a).subs(v, vt[0]))
print('hoop tension at terminal =', T, ' == mg(B0+Cz)/(pi a C)?', sp.simplify(T - m*g*(B0 + C*z)/(sp.pi*a*C)) == 0)
# half-ring force balance check numerically: force on half ring (phi in -pi/2..pi/2) x-component = 2T
half = sp.integrate((I*a*phat.cross(Bz*zhat))[0], (ph, -sp.pi/2, sp.pi/2))
print('half-ring radial resultant / 2 == T :', sp.simplify((half/2).subs(v, vt[0]) - T) == 0)
# numeric example
vals = {m: 2.0e-3, g: 9.8, R: 1.0e-3, a: 0.02, C: 5.0, B0: 0.4}
print('v_t numeric =', sp.N(vt[0].subs(vals)), 'm/s ; z where T=0 :', sp.N((-B0/C).subs(vals)))
print('T at z=0.10 m:', sp.N(T.subs(vals).subs(z, 0.10)), ' at z=-0.20:', sp.N(T.subs(vals).subs(z, -0.20)))
