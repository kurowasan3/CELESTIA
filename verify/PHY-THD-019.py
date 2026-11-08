# T20083-1 isothermal gas column: exact pressure difference, rigid body oscillation, flexible balloon
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

z, L, S, n, m, g, R, T, M, v, c0 = sp.symbols('z L S n m g R T M v c0', positive=True)
al = m*g/(R*T)
c = c0*sp.exp(-al*z)
# hydrostatic: dp/dz = -c m g with p = cRT
print('hydrostatic check', sp.simplify(sp.diff(c*R*T, z) + c*m*g))
c0sol = sp.solve(sp.Eq(S*sp.integrate(c, (z, 0, L)), n), c0)[0]
print('c(0) =', sp.simplify(c0sol), ' claim diff', sp.simplify(c0sol - n*al/(S*(1 - sp.exp(-al*L)))))
diff = sp.simplify((c.subs(z, 0) - c.subs(z, L)).subs(c0, c0sol))
print('c(0)-c(L) =', diff, ' claim nmg/(RTS) diff', sp.simplify(diff - n*m*g/(R*T*S)))
# equilibrium height
z0 = sp.log(m*v*c0/M)/al
print('z0 check', sp.simplify(c.subs(z, z0)*m*v - M))
# force near z0 : F = (c(z) m v - M) g ; derivative at z0
F = (c*m*v - M)*g
k_eff = sp.simplify(-sp.diff(F, z).subs(z, z0))
print('k_eff', k_eff, ' claim M g alpha diff', sp.simplify(k_eff - M*g*al))
vals = {m: 4.4e-2, g: 9.8, R: 8.3, T: 300}
alpha = float(al.subs(vals)); om = np.sqrt(9.8*alpha)
print('alpha', alpha, ' omega', om, ' period', 2*np.pi/om)
# numerical: nonlinear motion of rigid body in exponential atmosphere, small amplitude
Mv, vv, c0v = 1.0e-3, 1.0, 50.0  # kg, m^3, mol/m^3  -> rho0 = 2.2 kg/m^3
z0v = np.log(4.4e-2*vv*c0v/Mv)/alpha
f = lambda t, y: [y[1], (c0v*np.exp(-alpha*y[0])*4.4e-2*vv - Mv)*9.8/Mv]
s = solve_ivp(f, [0, 2000], [z0v + 1.0, 0], max_step=0.5, rtol=1e-10, atol=1e-10)
zz = s.y[0] - z0v
crossings = s.t[1:][(zz[:-1] < 0) & (zz[1:] >= 0)]
print('z0 [m]', z0v, ' simulated period', np.mean(np.diff(crossings)))
# flexible balloon: volume = n' R T / p(z); buoyancy = c(z) m g * volume = n' m g, independent of z
np_, p = sp.symbols("n' p", positive=True)
Vb = np_*R*T/(c*R*T)
print('balloon buoyancy', sp.simplify(c*m*g*Vb))
