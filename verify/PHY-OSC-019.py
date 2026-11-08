# P24-1 verification: plate A on vertical spring, block B on it; separation, reunion, inelastic collision
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

M, m, k, g, D, u, y = sp.symbols('M m k g D u y', positive=True)
d0 = (M + m)*g/k
# (2) energy: (1/2)(M+m)u^2 = (1/2)k D^2 - (M+m) g D   (spring natural length at y=0)
u2 = sp.solve(sp.Eq((M+m)*u**2/2, k*D**2/2 - (M+m)*g*D), u**2)[0]
print('u^2 =', sp.factor(u2))
# (4) u = pi g sqrt(M/k)
Dsol = sp.solve(sp.Eq(u2, sp.pi**2*g**2*M/k), D)
print('D solutions:', [sp.simplify(s) for s in Dsol])
Dexp = (M+m)*g/k*(1 + sp.sqrt(1 + sp.pi**2*M/(M+m)))
print('check D:', [sp.simplify(s - Dexp) for s in Dsol])
# arctan z < z  -> A returns first
z = sp.symbols('z', positive=True)
print('t_A = (2/w) atan(w u/g), t_B = 2u/g ; z - atan z derivative =', sp.simplify(sp.diff(z - sp.atan(z), z)))
# (5) collision
Vp = sp.simplify((M*u - m*u)/(M+m))
H = sp.simplify(Vp**2/(2*g)).subs(u, sp.pi*g*sp.sqrt(M/k))
print('H =', sp.factor(H))
heat = sp.simplify(M*u**2/2 + m*u**2/2 - (M+m)*Vp**2/2)
print('heat =', sp.factor(heat))

# ---------- numerical simulation with contact ----------
def simulate(Mv=2.0, mv=1.0, kv=50.0, gv=9.8, e_inelastic=True, tmax=6.0):
    Dv = float(Dexp.subs({M: Mv, m: mv, k: kv, g: gv}))
    w = np.sqrt(kv/Mv)
    # phase 1: together from y=-D at rest
    def rhs_tog(t, s): return [s[1], (-kv*s[0])/(Mv+mv) - gv]
    def contact_zero(t, s):  # contact force on B: N = m(g + a) where a = -k y/(M+m) - g -> N = -m k y/(M+m) -> zero at y=0
        return s[0]
    contact_zero.terminal = True; contact_zero.direction = 1
    s1 = solve_ivp(rhs_tog, [0, tmax], [-Dv, 0], events=contact_zero, rtol=1e-12, atol=1e-12)
    t_sep = s1.t_events[0][0]; useps = s1.y_events[0][0][1]
    print(f'separation at t={t_sep:.6f}, y={s1.y_events[0][0][0]:.2e}, u={useps:.6f}, expected u={np.pi*gv/w:.6f}')
    # phase 2: separate
    def rhs_sep(t, s):
        yA, vA, yB, vB = s
        return [vA, -kv*yA/Mv - gv, vB, -gv]
    def meet(t, s): return s[2] - s[0]
    meet.terminal = True; meet.direction = -1
    s2 = solve_ivp(rhs_sep, [t_sep, t_sep + tmax], [0, useps, 0, useps], events=meet, rtol=1e-12, atol=1e-12,
                   first_step=1e-6, max_step=1e-3)
    # skip trivial event at start: check time
    tm = s2.t_events[0]; ym = s2.y_events[0]
    idx = [i for i, tt in enumerate(tm) if tt - t_sep > 1e-4][0]
    tmeet = tm[idx]; st = ym[idx]
    print(f'meet after {tmeet - t_sep:.6f} (T_A={2*np.pi/w:.6f}) at yA={st[0]:.2e}, vA={st[1]:.5f}, vB={st[3]:.5f}')
    # check B always above A before meeting
    s2b = solve_ivp(rhs_sep, [t_sep, tmeet], [0, useps, 0, useps], rtol=1e-12, atol=1e-12, max_step=1e-3, dense_output=True)
    gap = s2b.y[2] - s2b.y[0]
    print('min gap before meeting (should be >=0):', gap[1:-1].min())
    # A returns to y=0 time (first downward crossing)
    ret = [s2b.t[i] for i in range(1, len(s2b.t)) if s2b.y[0][i-1] > 0 >= s2b.y[0][i]]
    print('A first returns at', ret[0] - t_sep, ' formula', 2/w*np.arctan(w*useps/gv), '; B returns at', 2*useps/gv)
    # phase 3: inelastic
    V = (Mv*st[1] + mv*st[3])/(Mv+mv)
    # continue together? contact force N = -m k y/(M+m) ... with y=0 and V>0 they separate immediately
    Hmax = V**2/(2*gv)
    print(f'V after collision {V:.6f}, B max height {Hmax:.6f}, formula {float(H.subs({M:Mv,m:mv,k:kv,g:gv})):.6f}')
    heatn = 0.5*Mv*st[1]**2 + 0.5*mv*st[3]**2 - 0.5*(Mv+mv)*V**2
    print('heat', heatn, 'formula', float(heat.subs({M:Mv,m:mv,k:kv,g:gv,u:useps})))
simulate()
simulate(Mv=3.0, mv=0.5, kv=80.0)
