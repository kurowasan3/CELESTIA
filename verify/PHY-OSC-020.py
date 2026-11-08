# P26-1 verification: ball in semicircular track on a free base, base latched when ball passes bottom
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

m, M, g, R, L = sp.symbols('m M g R L', positive=True)
# small oscillation with free base: relative coordinate s (horizontal offset of ball from bottom of track)
# ball: m(x2'' + s'') = -m g s/R ; base: M x2'' = + m g s/R
w1 = sp.sqrt(g/R*(1 + m/M))
T1 = 2*sp.pi/w1
print('T1 =', sp.simplify(T1))
# base displacement when s=0: CM fixed: m(x2 + s) + M x2 = m L -> x2 = m(L - s)/(m+M)
x2P = m*L/(m+M)
# relative speed at bottom from energy: (1/2) mu vrel^2 = m g L^2/(2R)
mu = m*M/(m+M)
vrel = sp.sqrt(m*g*L**2/R/mu)
vball = M/(m+M)*vrel
vbase = m/(m+M)*vrel
L2 = sp.simplify(vball*sp.sqrt(R/g))
print('ball speed at bottom:', sp.simplify(vball), ' new amplitude:', L2)
assert sp.simplify(L2**2 - L**2*M/(m+M)) == 0
lost = sp.simplify(M*vbase**2/2)
print('energy lost:', sp.factor(lost))
assert sp.simplify(lost - m**2*g*L**2/(2*(m+M)*R)) == 0
xmin = x2P - L2
print('leftmost x of ball:', sp.simplify(xmin))
ttot = T1/4 + sp.pi/2*sp.sqrt(R/g)
print('time to leftmost:', sp.simplify(ttot))

# ---- full nonlinear numerical simulation (angle theta, base x) ----
def sim(mv=1.0, Mv=3.0, gv=9.8, Rv=1.0, Lv=0.03):
    th0 = np.arcsin(Lv/Rv)
    # Lagrangian: ball pos = (X + R sin th, -R cos th)
    def rhs(t, y):
        X, th, Xd, thd = y
        # equations: (m+M) Xdd + m R cos th thdd - m R sin th thd^2 = 0
        #            m R cos th Xdd + m R^2 thdd + m g R sin th = 0
        A = np.array([[mv+Mv, mv*Rv*np.cos(th)], [mv*Rv*np.cos(th), mv*Rv**2]])
        b = np.array([mv*Rv*np.sin(th)*thd**2, -mv*gv*Rv*np.sin(th)])
        Xdd, thdd = np.linalg.solve(A, b)
        return [Xd, thd, Xdd, thdd]
    def bottom(t, y): return y[1]
    bottom.terminal = True; bottom.direction = -1
    s = solve_ivp(rhs, [0, 10], [0, th0, 0, 0], events=bottom, rtol=1e-12, atol=1e-14)
    tb = s.t_events[0][0]; X, th, Xd, thd = s.y_events[0][0]
    vb = Xd + Rv*np.cos(th)*thd
    # period check (full)
    s2 = solve_ivp(rhs, [0, 10], [0, th0, 0, 0], rtol=1e-12, atol=1e-14, dense_output=True, max_step=1e-3)
    # find successive maxima of theta via sign change of thd
    ts = s2.t; tds = s2.y[3]
    ups = [ts[i] for i in range(1, len(ts)) if tds[i-1] < 0 <= tds[i]]
    T1n = ups[1] - ups[0]
    T1f = 2*np.pi/np.sqrt(gv/Rv*(1+mv/Mv))
    print(f'T1 numeric {T1n:.6f}, formula {T1f:.6f}; t_bottom {tb:.6f} vs T1/4 {T1f/4:.6f}')
    print(f'base pos at bottom {X:.6e} formula {mv*Lv/(mv+Mv):.6e}; base speed {Xd:.6e}')
    # base max excursion
    print('base max X', s2.y[0].max(), ' 2mL/(m+M)=', 2*mv*Lv/(mv+Mv))
    # after latch: fixed base, ball pendulum with initial angular velocity
    thd_new = vb/(Rv*np.cos(th))
    def rhs_f(t, y): return [y[1], -gv/Rv*np.sin(y[0])]
    def turn(t, y): return y[1]
    turn.terminal = True; turn.direction = 1
    s3 = solve_ivp(rhs_f, [0, 10], [0.0, thd_new], events=turn, rtol=1e-12, atol=1e-14)
    thmin = s3.y_events[0][0][0]; t3 = s3.t_events[0][0]
    xball = X + Rv*np.sin(thmin)
    print(f'new amplitude (horiz) {-Rv*np.sin(thmin):.6e} formula {Lv*np.sqrt(Mv/(mv+Mv)):.6e}')
    print(f'leftmost x {xball:.6e} formula {mv*Lv/(mv+Mv) - Lv*np.sqrt(Mv/(mv+Mv)):.6e}')
    print(f'time to leftmost {tb + t3:.6f} formula {T1f/4 + np.pi/2*np.sqrt(Rv/gv):.6f}')
    # normal force check: ball stays on track (N>0) trivially for small amplitude
sim()
sim(mv=2.0, Mv=1.0)
