# M37-1 verification: two masses, three springs; middle spring length chosen so omega_anti = 2 omega_same
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

m, k, l, L, d, t = sp.symbols('m k l L d t', positive=True)
kc = k*l/L                     # spring constant inversely proportional to natural length (same material as outer)
K = sp.Matrix([[k + kc, -kc], [-kc, k + kc]])
ev = K.eigenvals()
print('eigen (k_eff):', ev)
w_same = sp.sqrt(k/m); w_anti = sp.sqrt((k + 2*kc)/m)
Lsol = sp.solve(sp.Eq(w_anti, 2*w_same), L)
print('(1) L =', Lsol, ' kc =', kc.subs(L, Lsol[0]))
w = w_same
xA = d/2*(sp.cos(w*t) + sp.cos(2*w*t)); xB = d/2*(sp.cos(w*t) - sp.cos(2*w*t))
kcv = kc.subs(L, Lsol[0])
# check equations of motion
eqA = sp.simplify(m*sp.diff(xA, t, 2) - (-k*xA - kcv*(xA - xB)))
eqB = sp.simplify(m*sp.diff(xB, t, 2) - (-k*xB - kcv*(xB - xA)))
print('(3) EOM residuals:', eqA, eqB, ' ICs:', xA.subs(t, 0), xB.subs(t, 0), sp.diff(xA, t).subs(t, 0), sp.diff(xB, t).subs(t, 0))
# (4) first time A at rest at x=0 : solve numerically over a fine grid
T1 = sp.pi*sp.sqrt(m/k)
print('(4) at t=pi sqrt(m/k): xA', sp.simplify(xA.subs(t, T1)), 'vA', sp.simplify(sp.diff(xA, t).subs(t, T1)),
      'xB', sp.simplify(xB.subs(t, T1)), 'vB', sp.simplify(sp.diff(xB, t).subs(t, T1)))
# middle-spring deformation and B max displacement
print('    xB - xA =', sp.simplify(sp.trigsimp(xB - xA)))
E0 = sp.Rational(1, 2)*k*d**2 + sp.Rational(1, 2)*kcv*d**2
print('(2) initial energy (work) =', sp.simplify(E0))

# numeric : check that no earlier instant has xA=0 and vA=0 simultaneously
mv, kv, dv = 1.0, 4.0, 0.1
kcv_ = 1.5*kv
rhs = lambda tt, y: [y[2], y[3], (-kv*y[0] - kcv_*(y[0]-y[1]))/mv, (-kv*y[1] - kcv_*(y[1]-y[0]))/mv]
ts = np.linspace(0, np.pi*np.sqrt(mv/kv)*1.0001, 200001)
s = solve_ivp(rhs, [0, ts[-1]], [dv, 0, 0, 0], t_eval=ts, rtol=1e-12, atol=1e-14)
score = np.abs(s.y[0])/dv + np.abs(s.y[2])/(dv*np.sqrt(kv/mv))
i = np.argmin(score[1000:]) + 1000
print('numeric: best (xA=0,vA=0) at t=%.6f (pi*sqrt(m/k)=%.6f) score %.2e; xB=%.6f vB=%.2e' % (ts[i], np.pi*np.sqrt(mv/kv), score[i], s.y[1][i], s.y[3][i]))
# list all near-zero vA instants with |xA| small before T1
cand = [(ts[j], s.y[0][j]) for j in range(1, len(ts)-1) if np.sign(s.y[2][j]) != np.sign(s.y[2][j+1])]
print('vA=0 instants (t, xA):', [(round(a, 5), round(b, 5)) for a, b in cand])
print('max xB before T1:', s.y[1].max(), ' (9d/16 =', 9*dv/16, ')')
