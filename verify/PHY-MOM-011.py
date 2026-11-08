# T20011-1 verification: pendulum hitting smooth wall; string (slack + jerk) vs rod (constraint impulse)
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

m, g, l, a, e, v0 = sp.symbols('m g l alpha e v0', positive=True)
v0e = sp.sqrt(2*g*l*sp.cos(a))
# B = (-l sin a, -l cos a); before: v0(-cos a, sin a); after (smooth wall, restitution e on x):
vB = sp.Matrix([e*v0*sp.cos(a), v0*sp.sin(a)])
u_out = sp.Matrix([-sp.sin(a), -sp.cos(a)])
print('(1) radial outward comp after =', sp.factor(sp.simplify(vB.dot(u_out))))
# (2) apex directly below O: horizontal travel l sin a during rise time v0 sin a/g
esol = sp.solve(sp.Eq(e*v0e*sp.cos(a)*(v0e*sp.sin(a)/g), l*sp.sin(a)), e)
print('(2) e =', [sp.simplify(s) for s in esol])
ee = esol[0]
# (3) at B' (l sin a, -l cos a): velocity (e v0 cos a, -v0 sin a); tangential unit (cos a, sin a)
Vt = sp.simplify((ee*v0*sp.cos(a)*sp.cos(a) - v0*sp.sin(a)**2))
print('(3) speed after jerk =', sp.simplify(Vt), ' == v0 cos2a/2 ?', sp.simplify(Vt - v0*sp.cos(2*a)/2) == 0)
Vr = sp.simplify(ee*v0*sp.cos(a)*sp.sin(a) + v0*sp.sin(a)*sp.cos(a))
print('    lost radial comp =', sp.simplify(Vr))
Hs = sp.simplify(-l*sp.cos(a) + (Vt**2/(2*g)).subs(v0, v0e))
print('    string: highest point rel. O =', sp.factor(Hs))
# (4) rod
Jr, Jw = sp.symbols('J_r J_w')
dp = m*(sp.Matrix([e*v0*sp.cos(a), -e*v0*sp.sin(a)]) - sp.Matrix([-v0*sp.cos(a), v0*sp.sin(a)]))
sol = sp.solve([sp.Eq(Jw - Jr*sp.sin(a), dp[0]), sp.Eq(-Jr*sp.cos(a), dp[1])], [Jr, Jw])
print('(4) J_r (outward push, >0 compress) =', sp.simplify(sol[Jr]), '  J_w =', sp.simplify(sol[Jw]))
Hr = sp.simplify((-l*sp.cos(a) + (ee*v0)**2/(2*g)).subs(v0, v0e))
print('(5) rod: highest point rel O =', sp.factor(Hr))
diff = sp.simplify(Hr - Hs)
print('    Hr - Hs =', sp.factor(diff))
for av in [10, 20, 30, 40, 44]:
    print('     a=%d: Hr-Hs = %.5f l' % (av, float(diff.subs({a: sp.rad(av), l: 1}))))
print('    alpha=30:', sp.nsimplify(Hs.subs({a: sp.pi/6, l: 1})), sp.nsimplify(Hr.subs({a: sp.pi/6, l: 1})))

# numeric simulation of string case, alpha=30deg, l=1, g=9.8
gv, lv, av = 9.8, 1.0, np.pi/6
ev = 1/(2*np.cos(av)**2)
v0v = np.sqrt(2*gv*lv*np.cos(av))
# circle motion phase from horizontal right (theta=pi/2 measured from downward vertical, ccw positive to the right?)
# use angle phi from downward vertical, position (l sin phi, -l cos phi). start phi=+pi/2, ends at phi=-alpha
s = solve_ivp(lambda t, y: [y[1], -gv/lv*np.sin(y[0])], [0, 5], [np.pi/2, 0], events=lambda t, y: y[0] + av, rtol=1e-11, atol=1e-12)
phi, phid = s.y_events[0][0]
vx, vy = lv*np.cos(phi)*phid, lv*np.sin(phi)*phid
print('numeric speed at B', np.hypot(vx, vy), 'theory', v0v)
x0, y0 = -lv*np.sin(av), -lv*np.cos(av)
vx, vy = -ev*vx, vy
def dist(t): return np.hypot(x0 + vx*t, y0 + vy*t - 0.5*gv*t*t) - lv
ts = np.linspace(1e-6, 2, 20000)
d = np.array([dist(t) for t in ts])
idx = np.where(d[1:] * d[:-1] < 0)[0][0]
tt = brentq(dist, ts[idx], ts[idx+1])
X, Y = x0 + vx*tt, y0 + vy*tt - 0.5*gv*tt*tt
VX, VY = vx, vy - gv*tt
print('retension at', X, Y, ' expected', lv*np.sin(av), -lv*np.cos(av))
print('max distance from O during flight <= l:', max(d[:idx+1]) <= 1e-9)
tang = np.array([-Y, X])/lv  # ccw tangent
Vt_num = VX*tang[0] + VY*tang[1]
print('tangential speed after jerk', Vt_num, ' theory', v0v*np.cos(2*av)/2)
