# M32-1 verification: A (frictionless, on spring) pushes B (kinetic friction mu) ; separation point shifted by friction
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

m, M, k, mu, g, d, x = sp.symbols('m M k mu g d x', positive=True)
a = (-k*x - mu*M*g)/(m+M)                 # common acceleration (moving right, +x)
N = sp.simplify(M*a + mu*M*g)            # B: M a = N - mu M g
print('(1) N =', sp.factor(N))
x1 = sp.solve(sp.Eq(N, 0), x)[0]
print('(2) x1 =', x1)
c = -mu*M*g/k; R = d - mu*M*g/k
print('    turning point together =', sp.simplify(c+R), ' need > x1 -> d >', sp.simplify(x1 + mu*M*g/k + mu*M*g/k))
v1sq = sp.simplify(k/(m+M)*(R**2 - (x1-c)**2))
print('(3) v1^2 =', sp.factor(v1sq))
# energy-method cross-check: (1/2)k d^2 - (1/2)k x1^2 - mu M g (d+x1) = (1/2)(m+M)v1^2
print('    energy check:', sp.simplify(sp.Rational(1,2)*k*d**2 - sp.Rational(1,2)*k*x1**2 - mu*M*g*(d+x1) - sp.Rational(1,2)*(m+M)*v1sq))
cos_t = sp.simplify(-(x1-c)/R)
print('    cos(omega t1) =', cos_t)
xB = x1 + v1sq/(2*mu*g)
xAmax = sp.sqrt(x1**2 + m*v1sq/k)
print('(4) xB_stop =', sp.simplify(xB))

# ---------- numerical ----------
def sim(mv, Mv, kv, muv, gv, dv):
    def together(t, y): return [y[1], (-kv*y[0] - muv*Mv*gv)/(mv+Mv)]
    def ev(t, y): return mv*(-kv*y[0]) / mv * 0 + (muv*mv*gv - kv*y[0])   # N*(m+M)/M
    ev.terminal = True; ev.direction = -1
    s = solve_ivp(together, [0, 100], [-dv, 0], events=ev, rtol=1e-11, atol=1e-13, max_step=1e-3)
    t1 = s.t_events[0][0]; xs, vs = s.y_events[0][0]
    # after separation
    tB = vs/(muv*gv)
    ts = np.linspace(0, 6*np.pi*np.sqrt(mv/kv) + tB, 200001)
    w = np.sqrt(kv/mv)
    xA = xs*np.cos(w*ts) + vs/w*np.sin(w*ts)
    xBt = np.where(ts < tB, xs + vs*ts - 0.5*muv*gv*ts**2, xs + vs**2/(2*muv*gv))
    gap = xBt - xA
    return t1, xs, vs, xA.max(), xBt[-1], gap[1:].min()

for vals in [(1.0, 2.0, 50.0, 0.3, 9.8, 0.5), (2.0, 1.0, 80.0, 0.2, 9.8, 0.3), (1.0, 1.0, 20.0, 0.1, 9.8, 0.4)]:
    mv, Mv, kv, muv, gv, dv = vals
    sub = {m: mv, M: Mv, k: kv, mu: muv, g: gv, d: dv}
    t1, xs, vs, xAm, xBs, gapmin = sim(*vals)
    t1f = float(sp.sqrt((m+M)/k).subs(sub)*(np.pi - np.arccos(float(-cos_t.subs(sub)))))
    print(vals, f'| x1 {xs:.6f}/{float(x1.subs(sub)):.6f} v1 {vs:.6f}/{float(sp.sqrt(v1sq).subs(sub)):.6f}'
          f' t1 {t1:.6f}/{t1f:.6f} xB {xBs:.6f}/{float(xB.subs(sub)):.6f} xAmax {xAm:.6f}/{float(xAmax.subs(sub)):.6f} min gap {gapmin:.3e}')
