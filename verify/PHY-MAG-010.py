# T20001-1 verification: charged particle in a rotating smooth tube in uniform B
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

m, q, B, w, r0, t = sp.symbols('m q B omega r0 t', positive=True)
wc = q*B/m
r = sp.Function('r')
# Tube rotates CLOCKWISE seen from above (theta_dot = -w), B along +z.
# Lab-frame polar eqs: m(r'' - r th'^2) = F_r ; m(r th'' + 2 r' th') = F_th
th_d = -w
v_r, v_th = sp.symbols('v_r v_th')
# Lorentz force q v x B, v = v_r rhat + v_th phihat ; rhat x zhat = -phihat ; phihat x zhat = rhat
F_r_L = q*v_th*B
F_th_L = -q*v_r*B
rr = sp.symbols('rr'); rd = sp.symbols('rd')
Fr = F_r_L.subs(v_th, rr*th_d)
radial = sp.simplify((Fr)/m + rr*th_d**2)   # r'' = ...
print('r_ddot =', sp.factor(radial))
Om2 = sp.factor(-radial/rr)
print('Omega^2 =', Om2, ' == w(wc-w)?', sp.simplify(Om2 - w*(wc-w)) == 0)
# tube normal force N (tangential): m*2*rd*th_d = F_th_L + N
N = sp.simplify(m*2*rd*th_d - F_th_L.subs(v_r, rd))
print('N =', sp.factor(N), ' == m(wc-2w) rd ?', sp.simplify(N - m*(wc-2*w)*rd) == 0)
# max of Omega
wmax = sp.solve(sp.diff(w*(wc-w), w), w)
print('omega at max Omega:', wmax, ' Omega_max =', sp.simplify(sp.sqrt((w*(wc-w)).subs(w, wmax[0]))))
# motor work from r0 (released at rest relative to tube) to r=0: dKE lab
Om = sp.sqrt(w*(wc-w))
KE0 = sp.Rational(1,2)*m*(w*r0)**2
KE1 = sp.Rational(1,2)*m*(Om*r0)**2
W = sp.factor(sp.simplify(KE1-KE0))
print('W(r0->0) =', W)
# check power integral: P = N * (r*w) * sign ; tube velocity at particle = -w r phihat; power of N = N*(-w r)
rt = r0*sp.cos(Om*t)
rdt = sp.diff(rt, t)
Nt = m*(wc-2*w)*rdt
P = Nt*(-w*rt)
Wint = sp.simplify(sp.integrate(P, (t, 0, sp.pi/(2*Om))))
print('W by integrating power =', sp.factor(Wint), ' match:', sp.simplify(Wint-W) == 0)

# numeric: free particle in B with initial velocity of tube point, omega = wc/2 -> stays on tube line
mv, qv, Bv, r0v = 1.0, 1.0, 2.0, 1.0
wcv = qv*Bv/mv
for wv in [wcv/2, 0.3*wcv]:
    def f(tt, y):
        x, yy, vx, vy = y
        ax = qv/mv*(vy*Bv); ay = qv/mv*(-vx*Bv)
        return [vx, vy, ax, ay]
    # start at (r0,0), velocity of tube point (clockwise): (0, -w r0)
    sol = solve_ivp(f, [0, 2*np.pi/wcv], [r0v, 0, 0, -wv*r0v], rtol=1e-10, atol=1e-12, dense_output=True)
    ts = np.linspace(0, 2*np.pi/wcv, 400)
    X, Y = sol.sol(ts)[0], sol.sol(ts)[1]
    # angle of tube = -w t ; perpendicular distance from tube line
    dperp = -X*np.sin(-wv*ts) + Y*np.cos(-wv*ts)
    print(f'w={wv/wcv:.2f} wc: max distance of free particle from tube line =', np.max(np.abs(dperp)))
    if abs(wv - wcv/2) < 1e-12:
        # circle of diameter r0 through C : center (r0/2,0)? check |pos - c| const
        R = np.hypot(X - X.mean(), Y - Y.mean())
        print('  circle radius spread', R.min(), R.max(), ' center', X.mean(), Y.mean())
        print('  lab speed spread', np.ptp(np.hypot(sol.sol(ts)[2], sol.sol(ts)[3])))

# numeric: constrained motion r'' = -Omega^2 r with confinement, check period
wv = 0.3*wcv
Omv = np.sqrt(wv*(wcv-wv))
sol = solve_ivp(lambda tt, y: [y[1], -qv*Bv*wv/mv*y[0] + wv**2*y[0]], [0, 2*np.pi/Omv], [r0v, 0], rtol=1e-11, atol=1e-12)
print('after one predicted period r =', sol.y[0, -1], ' rdot =', sol.y[1, -1])
