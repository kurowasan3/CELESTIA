# M36-1 verification: mass between two rubber cords; lower cord slackens; X-only centre lies outside its region
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

m, k, g, b = sp.symbols('m k g b', positive=True)
x = sp.symbols('x', real=True)        # displacement downward from equilibrium
sX = b + m*g/k                       # X stretch at equilibrium  (k sX = mg + k b)
print('(1) X stretch at eq:', sX, '; X taut at x=-sqrt3 b iff', sp.simplify(sX - sp.sqrt(3)*b), '> 0')
w2 = sp.sqrt(2*k/m); w1 = sp.sqrt(k/m)
A0 = sp.sqrt(3)*b
v = sp.sqrt(w2**2*(A0**2 - b**2))
print('(2) v at Y natural length:', sp.simplify(v))
# energy cross-check with method B (elastic + gravity energies), from x=-sqrt3 b to x=b
EX = lambda xx: sp.Rational(1, 2)*k*(sX + xx)**2
EY = lambda xx: sp.Rational(1, 2)*k*(b - xx)**2
Eg = lambda xx: -m*g*xx
lhs = EX(-A0) + EY(-A0) + Eg(-A0)
rhs = EX(b) + EY(b) + Eg(b) + sp.Rational(1, 2)*m*v**2
print('    energy method B check:', sp.simplify(lhs - rhs))
A1 = sp.sqrt((2*b)**2 + v**2/w1**2)
low = -b + A1
print('(3) lowest point below eq:', sp.simplify(low))
# cross-check lowest via method B : from x=-sqrt3 b (rest) to x=low (rest), Y slack at low
lhs2 = EX(-A0) + EY(-A0) + Eg(-A0)
xl = sp.symbols('xl', positive=True)
sol = sp.solve(sp.Eq(lhs2, EX(xl) + Eg(xl)), xl)
print('    method B lowest:', [sp.simplify(s) for s in sol])
tslack = 2*sp.acos(2*b/A1)/w1
print('(4) slack time:', sp.simplify(tslack))
tup = 2*(sp.pi - sp.acos(b/A0))/w2
T = sp.simplify(tslack + tup)
print('(5) period:', T)

# ---------- numerical ----------
def run(mv, kv, gv, bv):
    sXv = bv + mv*gv/kv
    def F(xx):
        f = mv*gv
        if sXv + xx > 0: f -= kv*(sXv + xx)
        if bv - xx > 0: f += kv*(bv - xx)
        return f
    rhs = lambda t, y: [y[1], F(y[0])/mv]
    ev_slack = lambda t, y: y[0] - bv
    ev_turn = lambda t, y: y[1]
    ev_turn.direction = 1   # velocity from + to - ... detect top turning (v from - to +) handled below
    s = solve_ivp(rhs, [0, 20*np.sqrt(mv/kv)], [-np.sqrt(3)*bv, 0.0], events=[ev_slack, ev_turn],
                  rtol=1e-11, atol=1e-13, max_step=1e-4)
    te = s.t_events[0]
    slack_time = te[1] - te[0]                # first interval x>b
    tops = s.t_events[1]                      # v crosses 0 upward -> top turning points
    period = tops[1] - tops[0] if len(tops) > 1 else np.nan
    return s.y[0].max(), slack_time, period, (sXv + s.y[0]).min()

for vals in [(1.0, 50.0, 9.8, 0.1), (0.3, 20.0, 9.8, 0.05)]:
    mv, kv, gv, bv = vals
    sub = {m: mv, k: kv, g: gv, b: bv}
    lowN, tsN, TN, minX = run(*vals)
    print(vals, f'| low {lowN:.6f}/{float(low.subs(sub)):.6f} slack {tsN:.6f}/{float(tslack.subs(sub)):.6f}'
          f' period {TN:.6f}/{float(T.subs(sub)):.6f} min X stretch {minX:.4f} (>0 ok)')
