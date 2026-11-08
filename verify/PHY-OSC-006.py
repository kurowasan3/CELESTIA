# M31-1 verification: pendulum in a cart whose acceleration switches from +a to -a
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

# ---------- symbolic ----------
al, g, l, m = sp.symbols('alpha g l m', positive=True)
gp = g/sp.cos(al)                      # effective gravity |g'|
# phase 1: start theta=0 at rest, equilibrium -alpha -> amplitude alpha -> extreme -2alpha
v1max = sp.sqrt(2*gp*l*(1-sp.cos(al)))
print('(1) vmax =', sp.simplify(v1max))
# phase 2: equilibrium +alpha, start -2alpha at rest -> amplitude 3alpha -> reaches +4alpha
Tmax = m*gp*(3-2*sp.cos(3*al))
print('(3) Tmax =', sp.simplify(Tmax))
# (4): switch at theta=-alpha with speed v1max, new equilibrium +alpha
cphi = sp.symbols('cphi')
sol = sp.solve(sp.Eq(sp.Rational(1, 2)*v1max**2 + gp*l*(1-sp.cos(2*al)), gp*l*(1-cphi)), cphi)[0]
print('(4) cos phi =', sp.simplify(sol), ' check:', sp.simplify(sol-(sp.cos(al)+sp.cos(2*al)-1)))

# ---------- numerical ----------
def run(alpha, switch, gv=9.8, lv=1.0):
    a = gv*np.tan(alpha)
    def rhs(A):
        return lambda t, y: [y[1], -(gv*np.sin(y[0]) + A*np.cos(y[0]))/lv]
    if switch == 'extreme':          # first time thetadot = 0 after start (theta = -2alpha)
        ev = lambda t, y: y[1]
        ev.direction = 1
    else:                            # first max |speed| : thetadot' = 0 -> theta = -alpha
        ev = lambda t, y: y[0] + alpha
        ev.direction = -1
    ev.terminal = True
    s1 = solve_ivp(rhs(a), [0, 50], [0.0, 0.0], events=ev, rtol=1e-11, atol=1e-12, max_step=1e-3)
    y0 = s1.y_events[0][0]
    s2 = solve_ivp(rhs(-a), [0, 20], y0, rtol=1e-11, atol=1e-12, max_step=1e-3)
    th, w = s2.y
    T = lv*w**2 + gv*np.cos(th) + a*np.sin(th)   # tension/m with A=-a : T/m = l w^2 + g cos - A sin
    return y0, th, T

for alpha in [0.2, 0.4, 0.5]:
    y0, th, T = run(alpha, 'extreme')
    print(f'alpha={alpha}: switch theta={y0[0]:.6f} (expect {-2*alpha:.6f}); max theta={th.max():.6f} (expect {4*alpha:.6f});'
          f' Tmax/m={T.max():.6f} (expect {9.8*(3-2*np.cos(3*alpha))/np.cos(alpha):.6f}); Tmin/m={T.min():.5f}')
    y0, th, T = run(alpha, 'center')
    phi_num = -(th.min()) + alpha     # angle from new effective vertical (+alpha) on the rear side
    print(f'   (4) switch theta={y0[0]:.6f}; cos phi num={np.cos(phi_num):.6f}, formula={np.cos(alpha)+np.cos(2*alpha)-1:.6f}; Tmin/m={T.min():.4f}')
# slack check: alpha slightly above pi/6 -> tension at start of phase 2 negative
alpha = np.pi/6 + 0.05
print('alpha>pi/6: initial tension/m in phase 2 =', 9.8/np.cos(alpha)*np.cos(3*alpha))
