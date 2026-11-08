# P19-1 verification: chord tunnel, SHM collisions at midpoint H
import sympy as sp
import numpy as np

m, M, e, g, R, L, v0 = sp.symbols('m M e g R L v0', positive=True)
w = sp.sqrt(g/R)

# (1) force along chord: -(m g / R) r cos(phi) = -(m g/R) x  -> SHM with omega = sqrt(g/R)
v0e = L*w
# (2) first collision: P (v0) hits Q (rest)
vP, vQ = sp.symbols('vP vQ')
sol = sp.solve([m*v0 - (m*vP + M*vQ), vP - vQ + e*(v0 - 0)], [vP, vQ], dict=True)[0]
vP1, vQ1 = sp.simplify(sol[vP]), sp.simplify(sol[vQ])
print('vP1 =', sp.factor(vP1), ' vQ1 =', sp.factor(vQ1))
assert sp.simplify(vP1 - (m - e*M)*v0/(m+M)) == 0
assert sp.simplify(vQ1 - (1+e)*m*v0/(m+M)) == 0
# Q amplitude = vQ1/omega <= L  with v0 = L w  -> (1+e) m <= m + M  -> M >= e m
print('cond:', sp.simplify((vQ1.subs(v0, v0e)/w)/L))

# (3) after T/2 both return to H with reversed velocities
uP, uQ = -vP1, -vQ1
sol2 = sp.solve([m*uP + M*uQ - (m*vP + M*vQ), vP - vQ + e*(uP - uQ)], [vP, vQ], dict=True)[0]
vP2, vQ2 = sp.factor(sol2[vP]), sp.factor(sol2[vQ])
print('vP2 =', vP2, ' vQ2 =', vQ2)
assert sp.simplify(vP2 + (m + e**2*M)*v0/(m+M)) == 0
assert sp.simplify(vQ2 + (1 - e**2)*m*v0/(m+M)) == 0

# (4) final: common velocity at H = +-m v0/(m+M), amplitude m L/(m+M); heat
Ei = m*v0e**2/2
Ef = (m+M)*(m*v0e/(m+M))**2/2
print('heat =', sp.factor(sp.simplify(Ei - Ef)))
assert sp.simplify(Ei - Ef - m*M*g*L**2/(2*(m+M)*R)) == 0

# ---- numerical event-driven simulation (exact SHM between collisions, check positions) ----
def simulate(mv=1.0, Mv=1.7, ev=0.6, gv=9.8, Rv=6.4e6, hv=3.0e6, ncol=60):
    Lv = np.sqrt(Rv**2 - hv**2); om = np.sqrt(gv/Rv)
    # numerically integrate with brute force small steps using full 2D gravity component along chord
    from scipy.integrate import solve_ivp
    def rhs(t, y):
        xP, vPn, xQ, vQn = y
        # acceleration along chord from uniform sphere: -(g/R)*r_vec . unit = -(g/R) x
        return [vPn, -(gv/Rv)*xP, vQn, -(gv/Rv)*xQ]
    def ev_col(t, y):
        return y[0] - y[2]
    ev_col.terminal = True; ev_col.direction = 0
    y = [-Lv, 0.0, 0.0, 0.0]; t = 0.0
    times, pos = [], []
    for k in range(ncol):
        s = solve_ivp(rhs, [t, t + 10*2*np.pi/om], y, events=ev_col, rtol=1e-11, atol=1e-6,
                      max_step=2*np.pi/om/400)
        # skip spurious event at start
        if len(s.t_events[0]) == 0:
            break
        tc = s.t_events[0][0]; yc = s.y_events[0][0]
        if tc - t < 1e-6:
            # start exactly at collision configuration (after collision positions equal): step a bit
            s = solve_ivp(rhs, [t, t + 1e-3*2*np.pi/om], y, rtol=1e-11, atol=1e-6)
            y = s.y[:, -1]; t = s.t[-1]
            continue
        xP, a, xQ, b = yc
        na = (mv*a + Mv*b - Mv*ev*(a - b))/(mv + Mv)
        nb = (mv*a + Mv*b + mv*ev*(a - b))/(mv + Mv)
        times.append(tc); pos.append(xP)
        y = [xP, na, xQ, nb]; t = tc
    return np.array(times), np.array(pos), Lv, om, y
tt, pp, Lv, om, yf = simulate()
T = 2*np.pi/om
print('first collision time / (T/4):', tt[0]/(T/4))
print('intervals / (T/2):', np.round(np.diff(tt[:8])/(T/2), 6))
print('collision positions / L:', np.round(pp[:8]/Lv, 8))
mv, Mv, ev = 1.0, 1.7, 0.6
vG = abs(mv*yf[1] + Mv*yf[3])/(mv+Mv)
print('final CM speed at H / (L*om):', vG/(Lv*om), ' expected m/(m+M)=', mv/(mv+Mv))
print('final relative speed / v0:', abs(yf[1]-yf[3])/(Lv*om))
# Q amplitude check after first collision
vQ1n = (1+ev)*mv/(mv+Mv)*Lv*om
print('Q amplitude / L =', vQ1n/om/Lv, '(<=1 since M>=e m)')
